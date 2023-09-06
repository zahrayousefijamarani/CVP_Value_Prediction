#include <inttypes.h>
#include "cvp.h"


#include "mypredictor.h"
#include <iostream>
#include <bitset>
#include <cassert>
#include <iterator>
#include <string>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>
#include <cstdlib>
using namespace std;

int seq_commit;


#define NOTLLCMISS (actual_latency < 150)
#define NOTL2MISS (actual_latency < 60)
#define NOTL1MISS (actual_latency < 12)
#define FASTINST (actual_latency ==1)
#define MFASTINST (actual_latency <3)

int MAX_CONF_LVP = 0;
int max_array = 2048;
long long int table[2048][3] ={{0}};
int k=0;
int l=0;

int all_preds = 0;
int correct_preds = 0;
int false_preds = 0;

long long int lookup (long long int address) {
	for (k=0; k< max_array; k++) { 
      if (table[k][0] == address) {
          if(table[k][2] >= MAX_CONF_LVP)
          {
                  return table[k][1];
                  
          }
          throw table[k][1];// show found but not used for prediction
      }
  }
	throw "error"; // show not found in table
}



void update (long long int address_here, long long int value_here, long long int pred_here) {
  int found =0;
  if(value_here != pred_here) {  
    for (k=0; k<max_array; k++) {
      if (table[k][0] == address_here) { // incorrect prediction
          table[k][1] = value_here;
          found=1;
				  if(table[k][2] != 0) 
            table[k][2] -= 1; // decrease confidence !!!!
      }
    }
    if (found==0) { // not found, add to table
      for (l=max_array-1; l>0; l--) {
        table[l][0] = table[l-1][0];
        table[l][1] = table[l-1][1];
			  table[l][2] = table[l-1][2]; 
      }
      table[0][0] = address_here;
      table[0][1] = value_here;
		  table[0][2] = 0;
    }
  }
	else // prediction was correct, increase confidence !!!!
	 {
		for (k=0; k<max_array; k++) {
      if (table[k][0] == address_here) { 
        table[k][2] += 1;
        return;
      }
    }
	}}


void getPredValue (ForUpdate * U, uint64_t & predicted_value)
{
	uint64_t pc = U->pc;
	long long int res = 0; 
	try{
		res = lookup(pc);  // found with high confidence
		predicted_value = res;
		U->pred = true;	
	}
	catch(long long int x){ } // found but confidence was low
	catch (...)  { } // not found

}


PredictionResult getPrediction(const PredictionRequest& req)
{
	PredictionResult result;

	ForUpdate *U;
	U = &Update[req.seq_no & (MAXINFLIGHT - 1)];
	U->pc = req.pc + req.piece;
	U->pred = false;

  getPredValue(U, result.predicted_value);
  result.speculate = (U->pred);
  return result;
}

void UpdateValuePred (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
  uint64_t pc = U->pc;
  long long int pred = 0;
  try{ 
    pred = lookup(pc); // predicted 
	  update(pc,actual_value,pred);
  }
  catch(long long int x){ update(pc,actual_value, x);  } // found but not predicted, last value is important(we save it in a global var)
  catch (...)  { update(pc,actual_value,actual_value-1);  } // not found not predicted, we should just add to table 'pred' is not important
}

void updatePredictor (uint64_t seq_no, uint64_t actual_addr, uint64_t actual_value, const mem_data_t & store_data, uint64_t actual_latency)
{
  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  if (U->todo == 1)
  {
    UpdateValuePred(U, actual_value, (int) actual_latency);
    U->todo = 0;
  }
  seq_commit = seq_no;
}

void speculativeUpdate (uint64_t seq_no,	// dynamic micro-instruction # (starts at 0 and increments indefinitely)
		   bool eligible,	// true: instruction is eligible for value prediction. false: not eligible.
		   uint8_t prediction_result,	// 0: incorrect, 1: correct, 2: unknown (not revealed)
		   // Note: can assemble local and global branch history using pc, next_pc, and insn.
		   uint64_t
		   pc, uint64_t next_pc, InstClass insn, uint8_t mem_size, bool is_pair, uint8_t piece,
		   // Note: up to 3 logical source register specifiers, up to 1 logical destination register specifier.
		   // 0xdeadbeef means that logical register does not exist.
		   // May use this information to reconstruct architectural register file state (using log. reg. and value at updatePredictor()).
		   uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst)
{

  // the framework does not really allow  to filter the predictions, so we predict every instruction
  ForUpdate *U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  
  if (eligible)
    {
      U->NbOperand =(src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef);
      U->todo = 1;
      U->INSTTYPE = insn;
      U->prediction_result = (prediction_result == 1);
    }

  bool isCondBr = insn == condBranchInstClass;
  bool isUnCondBr = insn == uncondIndirectBranchInstClass || insn == uncondDirectBranchInstClass;
  
  //path history 
  // just to have a longer history without (software) management
  if ((isCondBr) || (isUnCondBr))
  {
    if (pc != next_pc - 4)
    {
      for (int i = 7; i > 0; i--)
        gpath[i] = (gpath[i] << 1) ^ ((gpath[i - 1] >> 63) & 1);

      gpath[0] = (gpath[0] << 1) ^ (pc >> 2);
      gtargeth = (gtargeth << 1) ^ (next_pc >> 2);
    }
  }
}

void beginPredictor (int argc_other, char **argv_other)
{
}

void
endPredictor ()
{
	printf("---------------------------- zyj --------------------------\n");
	printf("All : %d\n", all_preds);
	printf("Correct : %d\n", correct_preds);
	printf("Incorrect : %d\n", false_preds);
}
