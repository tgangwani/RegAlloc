//===-- RegAllocBasic.cpp - Basic Register Allocator ----------------------===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
// This file defines the RABasic function pass, which provides a minimal
// implementation of the basic register allocator.
//
//===----------------------------------------------------------------------===//

#define DEBUG_TYPE "regalloc"
#include "llvm/CodeGen/Passes.h"
#include "AllocationOrder.h"
#include "LiveDebugVariables.h"
#include "RegAllocBase.h"
#include "Spiller.h"
#include "llvm/ADT/Statistic.h"
#include "llvm/CodeGen/LiveVariables.h"
#include "llvm/Analysis/AliasAnalysis.h"
#include "llvm/CodeGen/CalcSpillWeights.h"
#include "llvm/CodeGen/LiveIntervalAnalysis.h"
#include "llvm/CodeGen/LiveRangeEdit.h"
#include "llvm/CodeGen/LiveRegMatrix.h"
#include "llvm/CodeGen/LiveStackAnalysis.h"
#include "llvm/CodeGen/MachineBlockFrequencyInfo.h"
#include "llvm/CodeGen/MachineFrameInfo.h"
#include "llvm/CodeGen/MachineFunctionPass.h"
#include "llvm/CodeGen/MachineInstr.h"
#include "llvm/CodeGen/MachineLoopInfo.h"
#include "llvm/CodeGen/MachineRegisterInfo.h"
#include "llvm/CodeGen/RegAllocRegistry.h"
#include "llvm/CodeGen/SlotIndexes.h"
#include "llvm/CodeGen/VirtRegMap.h"
#include "llvm/PassAnalysisSupport.h"
#include "llvm/Support/Debug.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/Target/TargetMachine.h"
#include "llvm/Target/TargetRegisterInfo.h"
#include "llvm/Target/TargetInstrInfo.h"
#include <cstdlib>
#include <float.h>
#include <queue>
#include <iostream>
#include <map>

STATISTIC(NumStores, "Number of stores added");
STATISTIC(NumLoads , "Number of loads added");
STATISTIC(NumSpills, "Number of spilled live ranges");

// --helpers--
unsigned findMin(unsigned *array, unsigned total);

using namespace llvm;

static RegisterRegAlloc unitPRegAlloc("unitP", "unitProject register allocator",
                                      createUnitProjectRegisterAllocator);

static RegisterRegAlloc spillAllRegAlloc("spillAll", "spillAll register allocator",
                                      createSpillAllRegisterAllocator);

// Node in the interference graph
class RegisterLiveInterval
{
private:
  LiveInterval* LI;   // [A,B) 
  std::set<RegisterLiveInterval*> neighbors;  // set of neighbors
  std::set<unsigned> physRegInterferences;  // physical registers with which this node interferes
  unsigned virtReg; // the virtual register this node represents
  unsigned physReg; // the physical register this node is assigned to (color)
  bool isSpill;   // is this reg spilled to stack
  double spillWeight; // cost

public:
  RegisterLiveInterval(unsigned vreg);
  void setLiveInterval(LiveInterval* LI) { this->LI = LI; } 
  void markSpill(bool isSpill) { this->isSpill = isSpill; }
  bool isMarkedSpill() { return isSpill; }
  LiveInterval* getLiveInterval() { return LI; }
  void addNeighbor(RegisterLiveInterval* RLI) { neighbors.insert(RLI); }
  void deleteNeighbor(RegisterLiveInterval* RLI);
  std::set<RegisterLiveInterval*> getNeighbors() { return neighbors; }

  unsigned getVirtReg() { return virtReg; }
  unsigned getPhysReg() { return physReg; }
  void setPhysReg(unsigned preg)  { physReg = preg; }
  void setPhysRegInterferences(unsigned num)  { physRegInterferences.insert(num); }
  std::set<unsigned> getPhysRegInterferences()  { return physRegInterferences; }
  void setSpillWeight(double sw)  { spillWeight = sw; }
  double getSpillWeight() { return spillWeight; }
  bool isColored()  { return (physReg != ~0u); }
  void resetNode();
};

RegisterLiveInterval::RegisterLiveInterval(unsigned vreg): LI(NULL), virtReg(vreg), physReg(~0u), isSpill(false), spillWeight(DBL_MAX)
{}

// remove a neighbor
void RegisterLiveInterval::deleteNeighbor(RegisterLiveInterval* RLI)
{
  if(neighbors.count(RLI))
    neighbors.erase(RLI);
}

// clear out node attributes
void RegisterLiveInterval::resetNode()
{
  neighbors.clear();
  physRegInterferences.clear();
  physReg = ~0u;
  isSpill = false;
}

namespace {

class RAUnitP : public MachineFunctionPass
{
  // list of spill metrics evaluated
  enum spillMetric{ 
    costByDegree = 1, 
    costByDegreeSq, 
    cost,
    equalWeight,
    totalMetrics 
  };

private:
  bool spillAll;
  MachineFunction *MF;
  VirtRegMap *VRM;
  LiveRegMatrix *Matrix;
  MachineRegisterInfo *MRI;
  MachineLoopInfo *MLI;
  const TargetRegisterInfo *TRI;
  const TargetInstrInfo *TII;
  LiveIntervals *LIS;
  SlotIndexes *SI;
  RegisterClassInfo RegClassInfo;

  // live interval map for virtual registers
  // <virtual-reg number, Node-object>
  std::map<unsigned, RegisterLiveInterval*> iMap; 
  
  // set of physical registers with pre-existing LiveIntervals 
  std::set<unsigned> livePhysReg; 

  // set of virtual registers (nodes) which were spilled
  std::set<unsigned> spilledVirtReg;

  // StackSlotForVirtReg - Maps virtual regs to the frame index where these
  // values are spilled. Only used for Direct mapping
  IndexedMap<int, VirtReg2IndexFunctor> StackSlotForVirtReg;

public:
  RAUnitP(bool s);
  ~RAUnitP();

  virtual const char* getPassName() const {
    return spillAll ? "Spill All Register Allocator" : "Unit Project Register Allocator";
  }
  
  virtual void getAnalysisUsage(AnalysisUsage &AU) const;
  virtual bool runOnMachineFunction(MachineFunction &mf);
  
  // Returns the frame index 
  int getStackSpaceFor(unsigned VirtReg, const TargetRegisterClass *RC){
    int SS = StackSlotForVirtReg[VirtReg];
    if (SS != -1)
      return SS;          // Already has space allocated?

    // Allocate a new stack object for this spill location...
    int FrameIdx = MF->getFrameInfo()->CreateSpillStackObject(RC->getSize(),
                                                            RC->getAlignment());

    // Assign the slot.
    StackSlotForVirtReg[VirtReg] = FrameIdx;
    return FrameIdx;
  } 
  
  void buildRegisterLiveIntervals();
  void buildInterferenceGraph();
  void calculateSpillCosts();
  bool colorGraph();
  unsigned tryColorWithSpillMetric(unsigned sm);
  void insertSpills();
  void spillEverything();
  void handleCallSites();
  void rewriteRegisters();
  static char ID; 
};

char RAUnitP::ID = 0;

}  // end anon namespace

//constructor
RAUnitP::RAUnitP(bool s): MachineFunctionPass(ID), spillAll(s), StackSlotForVirtReg(-1) {
  initializeLiveDebugVariablesPass(*PassRegistry::getPassRegistry());
  initializeLiveIntervalsPass(*PassRegistry::getPassRegistry());
  initializeSlotIndexesPass(*PassRegistry::getPassRegistry());
  initializeRegisterCoalescerPass(*PassRegistry::getPassRegistry());
  initializeMachineSchedulerPass(*PassRegistry::getPassRegistry());
  initializeLiveStacksPass(*PassRegistry::getPassRegistry());
  initializeMachineDominatorTreePass(*PassRegistry::getPassRegistry());
  initializeMachineLoopInfoPass(*PassRegistry::getPassRegistry());
  initializeVirtRegMapPass(*PassRegistry::getPassRegistry());
  initializeLiveRegMatrixPass(*PassRegistry::getPassRegistry());
}

//destructor
RAUnitP::~RAUnitP(){
  
  // reclaim memory
  for(std::map<unsigned, RegisterLiveInterval*>::iterator I = iMap.begin(); I != iMap.end(); ++I)
    delete I->second;
}

//Every pass has a static char ID, a pointer to which is stored when
//addRequired/addPreserved is called
void RAUnitP::getAnalysisUsage(AnalysisUsage &AU) const {
  AU.setPreservesCFG();
  AU.addRequired<AliasAnalysis>();
  AU.addPreserved<AliasAnalysis>();
  AU.addRequired<LiveIntervals>();
  AU.addPreserved<LiveIntervals>();
  AU.addRequired<SlotIndexes>();
  AU.addPreserved<SlotIndexes>();
  AU.addRequired<LiveDebugVariables>();
  AU.addPreserved<LiveDebugVariables>();
  AU.addRequired<LiveStacks>();
  AU.addPreserved<LiveStacks>();
  AU.addRequired<MachineBlockFrequencyInfo>();
  AU.addPreserved<MachineBlockFrequencyInfo>();
  AU.addRequiredID(MachineDominatorsID);
  AU.addPreservedID(MachineDominatorsID);
  AU.addRequired<MachineLoopInfo>();
  AU.addPreserved<MachineLoopInfo>();
  AU.addRequired<VirtRegMap>();
  AU.addPreserved<VirtRegMap>();
  AU.addRequired<LiveRegMatrix>();
  AU.addPreserved<LiveRegMatrix>();  
  MachineFunctionPass::getAnalysisUsage(AU);
}

// Initialize the iMap structure
// Go over virtual registers in code and create RegisterLiveInterval objects
void RAUnitP::buildRegisterLiveIntervals()
{
  for(MachineFunction::iterator mfi = MF->begin(); mfi != MF->end(); mfi++)
    for(MachineBasicBlock::instr_iterator mbbii = mfi->instr_begin(); mbbii != mfi->instr_end(); mbbii++)
      for(unsigned i = 0, e = mbbii->getNumOperands(); i != e; ++i){ 
        MachineOperand &MO = mbbii->getOperand(i);
        
        // if operand not a register skip
        if(!MO.isReg()) continue;
        
        unsigned Reg = MO.getReg();
        
        // if operand is a physical register, add it to the set of physReg
        // having pre-existing live intervals. Register which alias are not
        // added multiple times
        if(TargetRegisterInfo::isPhysicalRegister(Reg))
        {
          bool aliasPresent = false;
          for(MCRegAliasIterator Units(Reg, TRI, true); Units.isValid(); ++Units)
            if(livePhysReg.count(*Units))
            {
              aliasPresent = true;
              break;
            }

         if(!aliasPresent)
          livePhysReg.insert(Reg);
        }

        // if operand is virtual register, create a node for it and add to iMap
        if(TargetRegisterInfo::isVirtualRegister(Reg) && !iMap.count(Reg)){
          RegisterLiveInterval* RLI = new RegisterLiveInterval(Reg);
          RLI->setLiveInterval(&LIS->getInterval(Reg));
          iMap.insert(std::pair<unsigned, RegisterLiveInterval*>(Reg, RLI));
        }
      }
}

// Create the interference graph between RegisterLiveInterval objects
void RAUnitP::buildInterferenceGraph()
{ 
  const TargetRegisterClass *RC;
  BitVector availableRegs;

  // calculate interference with physical registers. Also clear out virtual
  // register neighbors for a fresh graph built
  for(std::map<unsigned, RegisterLiveInterval*>::iterator I = iMap.begin(); I != iMap.end(); ++I)
  {
    I->second->resetNode();
    LiveInterval* VReg = I->second->getLiveInterval();
    RC = MRI->getRegClass(I->first);
    availableRegs = TRI->getAllocatableSet(*MF, RC);

    // use the livePhysReg set computed above to assign value to the
    // physRegInterferences attribute of the node. Aliasing is take care of here
    for(std::set<unsigned>::iterator II = livePhysReg.begin(); II != livePhysReg.end(); ++II)
      for(MCRegAliasIterator Units(*II, TRI, true); Units.isValid(); ++Units)
        if(availableRegs.test(*Units) && Matrix->checkRegUnitInterference(*VReg, *Units))
          I->second->setPhysRegInterferences(*Units);
  }

  // Naive O(n-square) approach to build interference between n nodes
  // Scope for improvement in future work
  for(std::map<unsigned, RegisterLiveInterval*>::iterator I = iMap.begin(); I != iMap.end(); ++I)
    for(std::map<unsigned, RegisterLiveInterval*>::iterator II = iMap.begin(); II != iMap.end(); ++II)
    {
     // don't want to capture interference with self 
     if(I->first == II->first) continue; 

     LiveInterval* G_LI = I->second->getLiveInterval();
     LiveInterval* E_LI = II->second->getLiveInterval();
     bool interference = false;
     
     // for all the definitions in LiveInterval-G, check if the any of the range
     // in the LiveInterval-E contains the definition
     for(LiveInterval::const_iterator J = G_LI->begin(); J != G_LI->end(); ++J)
     {
       for(LiveInterval::const_iterator JJ = E_LI->begin(); JJ != E_LI->end(); ++JJ)
       {
        SlotIndex G_start = J->start;
        if(JJ->contains(G_start)) 
        {
          interference = true;
          break;
        }
       }

       if(interference)
       {
         I->second->addNeighbor(II->second);
         II->second->addNeighbor(I->second);
         break;
       }
     }
    }
}

// Calculate spill costs for every node in the graph (i.e., every LiveInterval).
// Each def and use in the live range adds 10^(loopDepth of inst) to the cost
void RAUnitP::calculateSpillCosts()
{
	for(std::map<unsigned, RegisterLiveInterval*>::iterator it = iMap.begin(); it != iMap.end(); ++it)
	{
		double newSpillWeight = 0;
    unsigned vreg = it->first;
    
    // go over the def-use of virtual register
		for (MachineRegisterInfo::reg_iterator mriIter = MRI->reg_begin(vreg);
				MachineInstr *machInst = mriIter.skipInstruction(); )
		{
			unsigned loopDepth = MLI->getLoopDepth(machInst->getParent());
			
      if (loopDepth > 35) {
					loopDepth = 35; // Avoid overflowing the variable
			}
			
      std::pair<bool, bool> readWrite = machInst->readsWritesVirtualRegister(vreg);
			newSpillWeight += (readWrite.first + readWrite.second) * pow(10, loopDepth);
		}
		
    // update spillWeight attribute of node
    it->second->setSpillWeight(newSpillWeight);
	}
}

// Color the graph using a particular spillCost metric, passed as an input
// parameter. This function is called many time from ColorGraph()
unsigned RAUnitP::tryColorWithSpillMetric(unsigned sm)
{
  unsigned spillCount = 0;

  // Nodes is cloned from iMap. Nodes is modified throughout this procedure as
  // part of destroying and re-creating the interference graph.
  std::map<unsigned, RegisterLiveInterval*> Nodes = iMap;

  // the stack holds the virtual register number and RegisterLiveInterval
  // pointer in form of a pair
  std::stack<std::pair<unsigned, RegisterLiveInterval*> > nodeStack;
  
  const TargetRegisterClass *RC;
  BitVector availableRegs;

  // this loop goes on until all the nodes from the interference graph have been
  // inserted into the stack
  while(1)
  {
    bool Changed;

    // first pass - check which nodes can be trivially pushed to stack i.e.
    // degree is less than available colors
    do
    {
      Changed = false;
      for(std::map<unsigned, RegisterLiveInterval*>::iterator I = Nodes.begin(); I != Nodes.end();)
      {
        unsigned vreg = I->first;
        RegisterLiveInterval *RLI = I->second;
        
        // degree is calcuated by the total number of edges in interference
        // graph + the interferences with the physical register ranges. Note
        // that only those physical registers exist in getPhysRegInterferences()
        // which are in the allocatable set of this virtual register
        unsigned degree = RLI->getNeighbors().size() + RLI->getPhysRegInterferences().size();

        RC = MRI->getRegClass(vreg);
        availableRegs = TRI->getAllocatableSet(*MF, RC);
        std::set<unsigned> aliasedRegs;
        unsigned availableColors = 0;
        unsigned reg = availableRegs.find_first();

        // this 'for loop' calculates the number of colors available for the
        // node. Note that availableColors != availableRegs.count() in case the
        // node uses GPR8. e.g. if AL, AH are available to a node, these need to
        // be added as a single color since an interference with a node colored
        // with EAX renders both AL and AH unusable. Aliasing is taken care of here
        for(unsigned j = 0; j < availableRegs.count(); ++j)
        {
          bool overlap = false;
          for(std::set<unsigned>::iterator II = aliasedRegs.begin(); II != aliasedRegs.end(); ++II)
            if(TRI->regsOverlap(reg, *II))
            {
              overlap = true;
              break;
            }
        
          if(!overlap)
          {
            for (MCRegAliasIterator Units(reg, TRI, true); Units.isValid(); ++Units)
              aliasedRegs.insert(*Units);
            availableColors++;
          }
          reg = availableRegs.find_next(reg);
        }

        // is degree < K?
        if(degree < availableColors)
        {
          Changed = true;
          nodeStack.push(std::make_pair(vreg, RLI));

          // remove edges from the interference graph
          for(std::map<unsigned, RegisterLiveInterval*>::iterator II = Nodes.begin(); II != Nodes.end(); ++II)
            II->second->deleteNeighbor(RLI);

          // remove the node itself from the interference graph
          Nodes.erase(I++);
        }
        else
          ++I;
      }
    }while(Changed);

    if(Nodes.empty()) 
      break;

    // second pass - from the nodes that couldn't be trivially pushed to stack,
    // pick the one with least spillWeight
      
    // virtual register to erase
    unsigned eraseReg = 0;  
    double minWeight = DBL_MAX;

    // min spillWeight calculation
    for(std::map<unsigned, RegisterLiveInterval*>::iterator I = Nodes.begin(); I != Nodes.end(); ++I)
    {
      // LiveInterval not spillable, meaning that it was already spilled once
      // before
      if(!I->second->getLiveInterval()->isSpillable()) continue; 
      
      unsigned degree = I->second->getNeighbors().size() + I->second->getPhysRegInterferences().size();
      double spillWeight = I->second->getSpillWeight();

      switch(sm)
      {
        case RAUnitP::costByDegree :  {spillWeight = spillWeight / degree; break;}
        case RAUnitP::costByDegreeSq :  {spillWeight = spillWeight / (degree*degree); break;}
        case RAUnitP::cost :  break;
        case RAUnitP::equalWeight : {spillWeight = 0; break;}                      
        default : assert(false && "Unspecified metric used for spill-cost evaluation");                                            
      }
      
      if(spillWeight < minWeight)
      {
        minWeight = spillWeight;
        eraseReg = I->first;
      }
    }

    assert(eraseReg && "No node was selected for spilling");
    
    // remove edges from interference graph
    for(std::map<unsigned, RegisterLiveInterval*>::iterator I = Nodes.begin(); I != Nodes.end(); ++I)
      I->second->deleteNeighbor(Nodes[eraseReg]);
                              
    nodeStack.push(std::make_pair(eraseReg, Nodes[eraseReg]));
    Nodes.erase(eraseReg);
  }

  assert(Nodes.empty() && "We missed putting something into stack");

  // third pass - pop the nodes from stack and try to color
  while(!nodeStack.empty())
  {
    std::pair<unsigned, RegisterLiveInterval*> element = nodeStack.top();
    unsigned vreg = element.first; 
    
    // get a list of neighbors this node used to have in the interference graph
    std::set<RegisterLiveInterval*> neighbors = iMap[vreg]->getNeighbors();
    
    RC = MRI->getRegClass(vreg);
    availableRegs = TRI->getAllocatableSet(*MF, RC);
    AllocationOrder Order(vreg, *VRM, RegClassInfo);

    // for each neighbor, check if it has already been re-added from the stack
    // into the interference graph, and also if it has been colored
    for(std::set<RegisterLiveInterval*>::iterator I = neighbors.begin(); I != neighbors.end(); ++I)
    {
      unsigned v = (*I)->getVirtReg();
      
      // if the neighbor has a color, remove that from the set of colors
      // available to me. Aliasing is taken care of here
      if(Nodes.count(v) && Nodes[v]->isColored())
        for (MCRegAliasIterator Units(Nodes[v]->getPhysReg(), TRI, true); Units.isValid(); ++Units) 
          availableRegs.reset(*Units); 
    }

    if(!availableRegs.count())
    {
      iMap[vreg]->markSpill(true); 
      spillCount++;  
    }
    else{
      // if I have an available color, check if that color interferes with any of
      // the physical register interferences of this node
      unsigned physRegToAssign;

      do
      {
        physRegToAssign = Order.next();
      }while(!availableRegs.test(physRegToAssign));

      std::set<unsigned> physRegInterferences = iMap[vreg]->getPhysRegInterferences();
      bool overlap;

      do
      {
        overlap = false;
        
        // assign the physical register to this node. This can be overwritten if
        // there is a physical register interference later
        iMap[vreg]->setPhysReg(physRegToAssign);

        // iterate over the physical interferences of this node, and check if
        // the chosen color clashes with any. Aliasing is taken care by the
        // regsOverlap()
        for(std::set<unsigned>::iterator I = physRegInterferences.begin(); I != physRegInterferences.end(); ++I)
          if(TRI->regsOverlap(*I, physRegToAssign))
          {
            availableRegs.reset(physRegToAssign);
            overlap = true;
            break;
          }
        
        if(overlap)
        {
          if(!availableRegs.count())
          {
            iMap[vreg]->markSpill(true);
            spillCount++;
            iMap[vreg]->setPhysReg(~0u);
            break;
          } 

          // try out the next available color 
          do
          {
            physRegToAssign = Order.next();
          }while(!availableRegs.test(physRegToAssign));
        }
      }while(overlap);

    } 

    // pop the node from stack and re-insert to interference graph
    nodeStack.pop();
    Nodes.insert(std::pair<unsigned, RegisterLiveInterval*>(vreg, iMap[vreg]));
  }

  return spillCount;
}


// This function should check if the graph is K-colorable. If not, mark some
// RegisterLiveInterval as spillable 
bool RAUnitP::colorGraph()
{
  // spillCount stores the number of spills for each of the tried spillCost
  // metrics. findMin() is called later to return the least of these
  unsigned spillCount[RAUnitP::totalMetrics];

  for(unsigned I = 1; I < RAUnitP::totalMetrics; I++)
  {
    spillCount[I] = tryColorWithSpillMetric(I);
    if(spillCount[I] == 0)
      return true;
  
    buildInterferenceGraph();
  }
  
  unsigned min = findMin(spillCount, RAUnitP::totalMetrics);
  return (tryColorWithSpillMetric(min) == 0); 
}  
  
// -- This is our SPILLER Implementation --
// This function is responsible for changing the iMap structure according to the
// spill decisions. At the end of this function, iMap should contain "non-empty"
// & "relevant" Live Intervals. If a virtual register is split into (possibly many) 
// virtual registers, then the old register is removed from iMap and the new ones are 
// added with appropriate LiveIntervals
void RAUnitP::insertSpills()
{ 
  for(MachineFunction::iterator mfi = MF->begin(); mfi != MF->end(); mfi++){
    int num_spill = 0;
    for(MachineBasicBlock::instr_iterator mbbii = mfi->instr_begin(); mbbii != mfi->instr_end(); mbbii++){
      // we don't want to go over the spill instructions
      if (num_spill){
        num_spill--;
        continue;
      }
      for(unsigned i = 0, e = mbbii->getNumOperands(); i != e; ++i) {
        MachineOperand &MO = mbbii->getOperand(i);
        
        // if operand not a register skip
        if(!MO.isReg()) continue;
        
        unsigned Reg = MO.getReg();
        
        if(TargetRegisterInfo::isVirtualRegister(Reg)) 
        {
          RegisterLiveInterval* RLI = iMap[Reg];
          if(!RLI->isMarkedSpill()) continue;

          // since this register is being spilled, its LiveInterval won't be
          // required anymore. We would create new smaller LiveIntervals at the
          // def and use points of this register
          LIS->removeInterval(Reg);
        
          const TargetRegisterClass *RC = MRI->getRegClass(Reg);
          int FI = getStackSpaceFor(Reg, RC);

          // add this to the spilledVirtReg set for counting purposes only!
          spilledVirtReg.insert(Reg);

          if(MO.isDef())
          { 
            num_spill++;
          
            // properly placing the store-to-stack instruction
            MachineBasicBlock::instr_iterator MI = mbbii;
            for(int s = 0; s < num_spill; s++)
              MI++;

            int newReg = MRI->createVirtualRegister(RC);
            TII->storeRegToStackSlot(*mfi, MI, newReg, true, FI, RC, TRI); 
            MO.setReg(newReg);
            NumStores++;
          }  
          else if(MO.isUse())
          {
            int newReg = MRI->createVirtualRegister(RC);
            TII->loadRegFromStackSlot(*mfi, mbbii, newReg, FI, RC, TRI);
            MO.setReg(newReg);
            NumLoads++;
          } 
          else
            assert(true && "Operand not def or use");
        } 
      } 
    } 
    
    // repair SlotIndex after adding instructions to MF
    SI->repairIndexesInRange(mfi, mfi->begin(), mfi->end());
  } 

  // re-compute the Live-Intervals!
  // for any new virtual registers we may have created, add these to iMap
  for(unsigned a = 0, e = MRI->getNumVirtRegs(); a!=e; ++a) {
      unsigned reg = TargetRegisterInfo::index2VirtReg(a);

      // If virtReg is not marked to spill, don't touch
      if(iMap.count(reg) && !iMap[reg]->isMarkedSpill()) continue;

      // If virtReg is not live anywhere don't create the object
      if(LIS->createAndComputeVirtRegInterval(reg).empty()) 
      {
        // this is some old virtual register which got split into 
        // (possibly many) new registers. Delete it from iMap and reclaim memory
        if(iMap.count(reg))
        {
          delete iMap[reg];
          iMap.erase(reg);
        }
        LIS->removeInterval(reg);
        continue;
      }

      assert(!iMap.count(reg) && "An existing register should not reach here");

      RegisterLiveInterval* RLI = new RegisterLiveInterval(reg);
      iMap.insert(std::pair<unsigned, RegisterLiveInterval*>(reg, RLI));
      iMap[reg]->setLiveInterval(&LIS->getInterval(reg));
      iMap[reg]->getLiveInterval()->markNotSpillable();
  } 
}

// this function is responsible for saving the caller saved registers (EAX, EBX,
// ECX, EDX) before a CALL instruction and restoring them after the CALL
void RAUnitP::handleCallSites()
{
  if(spillAll) return;

  // create set of callee saved register
  std::set<unsigned> calleeSavedRegs;
  for (const uint16_t *I = TRI->getCalleeSavedRegs(MF); *I; ++I) 
    calleeSavedRegs.insert(*I);

  // iterate over all the CALL instructions in the program
  for(MachineFunction::iterator mfi = MF->begin(); mfi != MF->end(); mfi++){
    int save_restore = 0;
    for(MachineBasicBlock::instr_iterator mbbii = mfi->instr_begin(); mbbii != mfi->instr_end(); mbbii++)
    {
      if(save_restore){
        save_restore--;
        continue;
      }

      const MCInstrDesc &MCID = mbbii->getDesc();
      if (MCID.isCall())
      {
        // postCallInst is the MachineInst before which the load (restore) would
        // be inserted. Note that we don't handle the corner case when the CALL
        // instruction is BB terminator and postCallInst in invalid
        MachineBasicBlock::instr_iterator postCallInst = mbbii;
        postCallInst++;
        SlotIndex postCallInst_SI = SI->getInstructionIndex(postCallInst);
        SlotIndex callInst_SI = SI->getInstructionIndex(mbbii);
        
        for(std::map<unsigned, RegisterLiveInterval*>::iterator it = iMap.begin(); it != iMap.end(); ++it)
        {
          // the physical register assigned is callee saved, do nothing
          if(calleeSavedRegs.count(it->second->getPhysReg())) continue;

          // the register has to be live across the CALL instruction
          if(!(it->second->getLiveInterval()->liveAt(callInst_SI) && 
              it->second->getLiveInterval()->liveAt(postCallInst_SI))) continue;
          
          save_restore++;
          unsigned vreg = it->first;
          const TargetRegisterClass *RC = MRI->getRegClass(vreg);
          int FI = getStackSpaceFor(vreg, RC);
          
          // save-restore
          TII->storeRegToStackSlot(*mfi, mbbii, vreg, true, FI, RC, TRI); 
          TII->loadRegFromStackSlot(*mfi, postCallInst, vreg, FI, RC, TRI);
        }
      }

    }

    // repair SlotIndex after adding instructions to MF
    SI->repairIndexesInRange(mfi, mfi->begin(), mfi->end());
  }
}

// Graph is fully colored, assign Phys2Virt
void RAUnitP::rewriteRegisters()
{
  VRM->clearAllVirt();
  
  if(!spillAll)
    for(std::map<unsigned, RegisterLiveInterval*>::iterator it = iMap.begin(); it != iMap.end(); ++it)
    {
      unsigned vreg = it->first;
      unsigned preg = it->second->getPhysReg();
      assert(preg != ~0u && "We missed coloring some node properly");
      VRM->assignVirt2Phys(vreg, preg);
    }
  
  // code for Spill-All register allocator
  else
    for(std::map<unsigned, RegisterLiveInterval*>::iterator it = iMap.begin(); it != iMap.end(); ++it)
    {
      LiveInterval* VirtReg = it->second->getLiveInterval();
      unsigned vreg = it->first;
      AllocationOrder Order(vreg, *VRM, RegClassInfo);
      bool assigned = false;
      while (unsigned PhysReg = Order.next())
      {
        if(assigned) break;
      
        // Check for interference in PhysReg
        switch (Matrix->checkInterference(*VirtReg, PhysReg)) {
          case LiveRegMatrix::IK_Free:
            // PhysReg is available, allocate it.
            {
              // this updates VRM internally
              Matrix->assign(*VirtReg, PhysReg);
              assigned = true;
            }
          default: continue;
        }
      } 
    }
}

// this function marks all the LiveIntervals for spilling
void RAUnitP::spillEverything()
{
  for(std::map<unsigned, RegisterLiveInterval*>::iterator it = iMap.begin(); it != iMap.end(); ++it)
    it->second->markSpill(true);
}

bool RAUnitP::runOnMachineFunction(MachineFunction &mf) {
  DEBUG(dbgs() << "********** UNIT PROJECT REGISTER ALLOCATION **********\n"
               << "********** Function: "
               << mf.getName() << '\n');

  MF = &mf;
  VRM = &getAnalysis<VirtRegMap>();
  Matrix = &getAnalysis<LiveRegMatrix>();
  MRI = &MF->getRegInfo();
  TRI = MF->getTarget().getRegisterInfo();
  TII = MF->getTarget().getInstrInfo();
  LIS = &getAnalysis<LiveIntervals>();
  MLI = &getAnalysis<MachineLoopInfo>();
  SI = &getAnalysis<SlotIndexes>();
  MRI->freezeReservedRegs(mf);
  RegClassInfo.runOnMachineFunction(mf);
  
  iMap.clear();
  StackSlotForVirtReg.resize(MRI->getNumVirtRegs());
  
  buildRegisterLiveIntervals();
  
  if(spillAll)
  {
    spillEverything();
    insertSpills();
  }
  else
    while(1)
    { 
      buildInterferenceGraph();
      calculateSpillCosts();
      if(colorGraph())
        break;
      insertSpills();
    }
  
  handleCallSites();
  rewriteRegisters();

  // report statistics
  NumSpills = spilledVirtReg.size();

  StackSlotForVirtReg.clear();
  return true;
}

FunctionPass* llvm::createUnitProjectRegisterAllocator()
{
  return new RAUnitP(false);
}

FunctionPass* llvm::createSpillAllRegisterAllocator()
{
  return new RAUnitP(true);
}

// -- helpers --
unsigned findMin(unsigned *array, unsigned total)
{
 unsigned min = ~0u;
 unsigned min_pos = ~0u;

 for(unsigned I = 1; I < total; I++)
   if(array[I] < min)
   {
    min = array[I];
    min_pos = I;
   }

 return min_pos;
}
