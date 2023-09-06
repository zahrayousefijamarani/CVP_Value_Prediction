/*
 * MIT License
 * 
 * Copyright (c) 2019 Chirag Sakhuja
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
// 4.030
#define STRIDEON 
#define VTAGEON	
#define HCVPON	
#define VSEPON

int NSTRIDE, NHCVP, NVTAGE, NVSEP, NTOTAL, NPRED;
int NPREDSTRIDE, NPREDHCVP, NPREDVTAGE, NPREDVSEP; 

int NLOADSTRIDE, NLOADHCVP, NLOADVTAGE, NLOADVSEP, NLOADTOTAL, NLOADPRED;
int NLOADPREDSTRIDE, NLOADPREDHCVP, NLOADPREDVTAGE, NLOADPREDVSEP; 

#include "mypredictor.h"

/*
 * Much of this code is derived from the EVES CVP submission with the following modifications:
 * Code for VSEP was added
 * The code added for HCVP is clearly delimited with comments. Most of the logic to implement
 * HCVP is in mypredictor.h.
 */


#define dfcm

#ifdef dfcmplus
#include "dfcmplus.cc"
#endif


int seq_commit;

#define NOTLLCMISS (actual_latency < 150)
#define NOTL2MISS (actual_latency < 60)
#define NOTL1MISS (actual_latency < 12)
#define FASTINST (actual_latency ==1)
#define MFASTINST (actual_latency <3)

/** BEGIN FOR HCVP **/
std::unordered_map < uint64_t, InstInfo > seq_no_to_info;
ValueCorrelationPredictor
  VP;
/** END FOR HCVP **/

void
getPredStride (ForUpdate * U, uint64_t & predicted_value, uint64_t seq_no);
bool
strideupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency,
		  int stride);
bool
StrideAllocateOrNot (ForUpdate * U, uint64_t actual_value,
		     int actual_latency);
void
UpdateStridePred (ForUpdate * U, uint64_t actual_value, int actual_latency);
void
getPredVtage (ForUpdate * U, uint64_t & predicted_value);
bool
vtageupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency);
void
UpdateVtagePred (ForUpdate * U, uint64_t actual_value, int actual_latency);

void
getPredVSEP (ForUpdate * U, uint64_t & predicted_value, uint64_t seq_no);
bool
VSEPupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency);
void
UpdateVSEPPred (ForUpdate * U, uint64_t actual_value, int actual_latency);

// ---------------- fcm --------------------
#ifdef dfcm

static DFCM_Predictor predictor;

// Global branch and path history
static uint64_t ghr = 0, phist = 0;

// Load/Store address history
static uint64_t addrHist = 0;

int fcm_conf=9;//6

// the gash function is the Shift-Xor hash function
uint64_t gash(uint64_t a , uint64_t b, uint64_t c, uint64_t d)
{
	uint64_t output = 0;
	output ^= a;
	output ^= (b>>1);
	output ^= (c>>2);
	output ^= (d>>3);
	return output;
	
}


bool
getPredDFCM (ForUpdate * U, uint64_t & predicted_value, const PredictionRequest& req)
{
  // Calculating the index.
  // Instructions are 4-bytes, so shift PC by 2.
  uint64_t index = ((req.pc >> 2) ^ req.piece) & predictor.indexMask;



  // Accessing the lastValue, strideList, tag from the Level_1_Table.
  deque<uint64_t> strideList = predictor.Level_1_Table[index].strideList;
  uint64_t lastValue = predictor.Level_1_Table[index].lastValue;
  
  
  uint64_t a=strideList[0] & predictor.indexMask;
  uint64_t b=strideList[1] & predictor.indexMask;
  uint64_t c=strideList[2] & predictor.indexMask;
  uint64_t d=strideList[3] & predictor.indexMask;
  //uint64_t e=strideList[4] & predictor.indexMask;
  
  // Calculating the stride_index using the gash function
  uint64_t stride_index = gash(a,b,c,d) & predictor.indexMask;


  //Accessing the stride from the stride_Table
  uint64_t stride = predictor.Level_2_Table[stride_index].stride;

  predicted_value = lastValue + stride;
  uint8_t use_pred = (predictor.Level_2_Table[stride_index].conf >=fcm_conf);
  
  // Speculate using the prediction only if confidence is high enough
  return use_pred;

}

void
UpdateDFCM (ForUpdate * U, uint64_t actual_value, int actual_latency, uint64_t seq_no)
{

  std::deque<DFCM_Predictor::InflightInfo> & inflight = predictor.inflightPreds;

  // If we have value predictions waiting for corresponding update
  if(inflight.size() && seq_no == inflight.back().seqNum)
  {     
   
	uint64_t new_stride =   actual_value - predictor.Level_1_Table[inflight.back().index].lastValue; 
	
    // If there are not other predictions corresponding to this PC in flight in the pipeline, 
    // we make sure the DFCM stride history used to predict is clean by copying the architectural value
    // history in it.

    if(--predictor.Level_1_Table[inflight.back().index].inflight == 0)
    {
		predictor.Level_1_Table[inflight.back().index].strideList.pop_back();	
		predictor.Level_1_Table[inflight.back().index].strideList.push_front(predictor.Level_1_Table[inflight.back().index].lastValue);	
		predictor.Level_1_Table[inflight.back().index].lastValue = actual_value;
    }

	// If the predicted stride is equal to new stride , increase confidence.
    if(predictor.Level_2_Table[inflight.back().stride_Index].stride == new_stride)
	{
	     predictor.Level_2_Table[inflight.back().stride_Index].conf = 
	     std::min(predictor.Level_2_Table[inflight.back().stride_Index].conf + 1,8);
    }

	// if confidence is zero ,change the stride value in the predictor to new_stride.
    else if(predictor.Level_2_Table[inflight.back().stride_Index].conf == 0)
      predictor.Level_2_Table[inflight.back().stride_Index].stride = new_stride;

	// If prediction was wrong ,but confidence is not zero , then decrease confidence.
    else
    { 
	predictor.Level_2_Table[inflight.back().stride_Index].conf--;
        if (predictor.Level_2_Table[inflight.back().stride_Index].conf > 0)
	  predictor.Level_2_Table[inflight.back().stride_Index].conf -=1;
	if (predictor.Level_2_Table[inflight.back().stride_Index].conf > 0)
          predictor.Level_2_Table[inflight.back().stride_Index].conf -=1;
//        if (predictor.Level_2_Table[inflight.back().stride_Index].conf > 0 && U->INSTTYPE == loadInstClass)
//          predictor.Level_2_Table[inflight.back().stride_Index].conf -=1;

    }		
    inflight.pop_back();
  }


}
#endif
// ----------------------------------------



// --------------- lvp ------------------------------

int MAX_CONF_LVP = 6;
int max_array = 2048;
long long int table[2048][3] ={{0}};
int k=0;
int l=0;

long long int lookup (long long int address) {
        for (k=0; k< max_array; k++) {
                if (table[k][0] == address) {
                        if(table[k][2] >= MAX_CONF_LVP)
                        {
                               return table[k][1];

                        }
                      throw table[k][1];// show found but not used for prediction^M
                }
        }
        throw "error"; // show not found in table^M
}

void update (long long int address_here, long long int value_here, long long int pred_here) {
         int found =0;
         if(value_here != pred_here) {
                for (k=0; k<max_array; k++) {
                        if (table[k][0] == address_here) { // incorrect prediction^M
                                table[k][1] = value_here;
                                found=1;
                                table[k][2] -= 2;
                                if(table[k][2] < 0) { // decrease confidence !!!!
                                        table[k][2] = 0;
                                }
                        }
                }
                if (found==0) { // not found, add to table^M
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
         else // prediction was correct, increase confidence !!!!^M
         {
                for (k=0; k<max_array; k++) {
                        if (table[k][0] == address_here) {
                                table[k][2] += 1;
                                if(table[k][2] >= (MAX_CONF_LVP * 2-1)) table[k][2] = 2 * MAX_CONF_LVP-1;
                                return;
                        }
                }

        }}



// -------------- lvp -----------------------------------





PredictionResult getPrediction(const PredictionRequest& req)
{
  PredictionResult result;

    /** BEGIN FOR HCVP **/
  InstInfo
    inst_info
  {
    req.seq_no,
      req.pc,
  req.piece};
  seq_no_to_info[req.seq_no] = inst_info;
    /** END FOR HCVP **/

  ForUpdate *
    U;
  U = &Update[req.seq_no & (MAXINFLIGHT - 1)];
  U->pc = req.pc + req.piece;
  U->predstride = false;
  U->predhcvp = false;
  U->predvtage = false;
  U->predvsep = false;
   
#ifdef  HCVPON
    /** BEGIN FOR HCVP **/
  // A. Seznec: In practice better to give priority to HCVP on E-stride
  {

    auto
      corr_prediction = VP.getPrediction (req.seq_no, inst_info);
    if (std::get < 0 > (corr_prediction))
      {
	result.predicted_value = std::get < 1 > (corr_prediction);
	U->predicted_value = result.predicted_value;
	U->predhcvp = true;
      }
  }
      /** END FOR HCVP **/
#endif



#ifdef STRIDEON
  getPredStride (U, result.predicted_value, req.seq_no);
#endif

#ifdef VTAGEON
  getPredVtage (U, result.predicted_value);
#endif


//#ifdef VSEPON
//  getPredVSEP (U, result.predicted_value, req.seq_no);
//#endif

U->predicted_value = result.predicted_value;

  U->predicted = (U->predstride || U->predhcvp || U->predvtage || U->predvsep);
  result.speculate = U->predicted;

  if(result.speculate == false) // && (U->INSTTYPE == loadInstClass))  
{

#ifdef dfcm 
  	result.speculate = getPredDFCM(U, result.predicted_value, req);
#endif
        try{
        long long int res = lookup(U->pc);  // found with high confidence^M
        result.predicted_value = res;
        result.speculate = true;
        }
        catch(long long int x){ } // found but confidence was low^M
        catch (...)  { } // not found^M
  }

#ifdef VSEPON
  getPredVSEP (U, result.predicted_value, req.seq_no);
#endif

  return result;
}

/////////Update of  VTAGE
// function determining whether to  update or not confidence on a correct prediction

void
updatePredictor (uint64_t
                 seq_no,
                 uint64_t
                 actual_addr, uint64_t actual_value, const mem_data_t & store_data, uint64_t actual_latency)
{
  (void) actual_addr;

  ForUpdate *
    U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];

#ifdef HCVPON
    /** BEGIN FOR HCVP **/
  auto
    seq_no_to_info_search = seq_no_to_info.find (seq_no);

  if (seq_no_to_info_search->second.eligible)
    {
//         VP.updatePredictor(seq_no, actual_value);
// A. Seznec
#define UPDATECHCVP ((random () & (((1 << (2*FASTINST+NOTL2MISS+NOTL1MISS)- 1))))==0)
      bool
	UPDATECONFHCVP;
      switch (U->INSTTYPE)
	{
	case aluInstClass:

	case fpInstClass:
	case slowAluInstClass:
	case undefInstClass:
	case storeInstClass:
	  UPDATECONFHCVP = UPDATECHCVP;
	  break;
	case loadInstClass:
	  UPDATECONFHCVP = UPDATECHCVP || UPDATECHCVP;
	  break;
	case uncondIndirectBranchInstClass:
	  UPDATECONFHCVP = true;
	  break;
	default:
	  UPDATECONFHCVP = false;
	};
      VP.updatePredictor (seq_no, actual_value, UPDATECONFHCVP);
    }
  seq_no_to_info.erase (seq_no_to_info_search);
#endif
    /** END FOR HCVP **/
if (U->todo)
    {
         //        if ((random() & 7)==0)
         {
              if (!U->predvtage) if (!U->predvsep) if(U->predhcvp)
                                              {
                                                   if (actual_value != U->predicted_value) {if (ThresHCVP<32) ThresHCVP++;} else if ((random() & 127)==0)if (ThresHCVP> 14) ThresHCVP--;
                                              }
  if (U->predvtage)      if (actual_value != U->predicted_value) {if (ThresVTAGE<32) ThresVTAGE++;} else if ((random() & 511)==0)if (ThresVTAGE> 14) ThresVTAGE--;
if (!U->predvtage)  if (U->predvsep)      if (actual_value != U->predicted_value) {if (ThresVSEP<32) ThresVSEP++;} else if ((random() & 511)==0)if (ThresVSEP> 14) ThresVSEP--; 
         }
         
  if (U->predvtage)
	{
	  NVTAGE += (U->predvtage) & (actual_value != U->predicted_value);
	  NPREDVTAGE++;
	}
  else
       if (U->predvsep)
	{

	  NVSEP += (U->predvsep) & (actual_value != U->predicted_value);
	  NPREDVSEP++;
	} else
            
           if (U->predhcvp)
	{
	  NHCVP += (U->predhcvp) && (actual_value != U->predicted_value);
	  NPREDHCVP++;
	}
           else
                if (U->predstride)
	{
	  NSTRIDE += (U->predstride) & (actual_value != U->predicted_value);
	  NPREDSTRIDE++;
	}
if (U->INSTTYPE== loadInstClass){if (U->predvtage)
	{
	  NLOADVTAGE += (U->predvtage) & (actual_value != U->predicted_value);
	  NLOADPREDVTAGE++;
	}
      else
       if (U->predvsep)
	{

	  NLOADVSEP += (U->predvsep) & (actual_value != U->predicted_value);
	  NLOADPREDVSEP++;
	} else
            
           if (U->predhcvp)
	{
	  NLOADHCVP += (U->predhcvp) && (actual_value != U->predicted_value);
	  NLOADPREDHCVP++;
	}
           else
                if (U->predstride)
	{
	  NLOADSTRIDE += (U->predstride) & (actual_value != U->predicted_value);
	  NLOADPREDSTRIDE++;
	}
  }
  
               
    }
  if (U->todo == 1)
    {
      //just to force allocations and update on both predictors
      U->prediction_result = 0;
#ifdef STRIDEON
      UpdateStridePred (U, actual_value, (int) actual_latency);
#endif
#ifdef VSEPON
      UpdateVSEPPred (U, actual_value, (int) actual_latency);
#endif
#ifdef VTAGEON
      UpdateVtagePred (U, actual_value, (int) actual_latency);
#endif


#ifdef dfcm
      UpdateDFCM(U, actual_value, (int) actual_latency, seq_no);
#endif

      try{
        long long int pred = lookup(U->pc); // predicted ^M
        update(U->pc,actual_value,pred);
      }
      catch(long long int x){ update(U->pc,actual_value, x);  } // found but not predicted, last value is important(we save it in a global var)^M
      catch (...)  { update(U->pc,actual_value,actual_value-1);  } // not found not predicted, we should just add to table 'pred' is not important^M
//      }

      U->todo = 0;
    }
  seq_commit = seq_no;
}

void
speculativeUpdate (uint64_t seq_no,     // dynamic micro-instruction # (starts at 0 and increments indefinitely)
                   bool eligible,       // true: instruction is eligible for value prediction. false: not eligible.
                   uint8_t prediction_result,   // 0: incorrect, 1: correct, 2: unknown (not revealed)
                   // Note: can assemble local and global branch history using pc, next_pc, and insn.
                   uint64_t
                   pc, uint64_t next_pc, InstClass insn, uint8_t mem_size, bool is_pair, uint8_t piece,
                   // Note: up to 3 logical source register specifiers, up to 1 logical destination register specifier.
                   // 0xdeadbeef means that logical register does not exist.
                   // May use this information to reconstruct architectural register file state (using log. reg. and value at updatePredictor()).
                   uint64_t src1, uint64_t src2, uint64_t src3, uint64_t dst)
{

// ---------------------------- fcm ----------
#ifdef dfcm

  uint64_t index = ((pc >> 2) ^ piece) & predictor.indexMask;
  deque<uint64_t> strideList = predictor.Level_1_Table[index].strideList;


  uint64_t a=strideList[0] & predictor.indexMask;
  uint64_t b=strideList[1] & predictor.indexMask;
  uint64_t c=strideList[2] & predictor.indexMask;
  uint64_t d=strideList[3] & predictor.indexMask;
 // uint64_t e=strideList[4] & predictor.indexMask;
  
  // Calculating the stride_index using the gash function
  uint64_t stride_index = gash(a,b,c, d) & predictor.indexMask;


// if(isPredictable)
  if(eligible)
  {

	 // pusihng the instruction's inflightInfo into InflightPreds.
	 // increasing the inflight number for the given PC to indicate 
	 // one more instruction inflight which is yet to be committed.
	predictor.inflightPreds.push_front({seq_no, index, stride_index});
	predictor.Level_1_Table[index].inflight++;
    
//	return;
  }
  else { 

  
  // At this point, any branch-related information is architectural, i.e.,
  // updating the GHR/LHRs is safe.
  bool isCondBr = insn == condBranchInstClass;
  bool isIndBr = insn == uncondIndirectBranchInstClass;
  
  // Infrastructure provides perfect branch prediction.
  // As a result, the branch-related histories can be updated now.
  if(isCondBr)
    ghr = (ghr << 1) | (pc != next_pc + 4);

  if(isIndBr)
	phist = (phist << 4) | (next_pc & 0x3);
  }
#endif
//}
// -------------------------------------------


  (void) piece;
  (void) dst;
 bool
    isCondBr = insn == condBranchInstClass;
  bool
    isUnCondBr = insn == uncondIndirectBranchInstClass
    || insn == uncondDirectBranchInstClass;
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
// the framework does not really allow  to filter the predictions, so we predict every instruction    /** BEGIN FOR HCVP **/

  auto
    seq_no_to_info_search = seq_no_to_info.find (seq_no);
  seq_no_to_info_search->second.type = insn;
  seq_no_to_info_search->second.eligible = eligible;

  if ((insn == condBranchInstClass) || (insn == uncondIndirectBranchInstClass)
      || (insn == uncondDirectBranchInstClass))
    {
//A. Seznec
// I did modify a little bit to use path history instead of branch history         
         uint64_t
	outcome = 0;
      if ((pc + 4) != next_pc)
	{
             outcome = (pc ^ (pc >>2) ^ (pc >>6));
          VP.shiftBHR (outcome);
	}

    }
    /** END FOR HCVP **/

  ForUpdate *
    U;
  U = &Update[seq_no & (MAXINFLIGHT - 1)];
  
  LastMispVT++;
  LastMispHCVP++;
  LastMispVSEP++;

  if (eligible)
    {
      U->NbOperand =
	(src1 != 0xdeadbeef) + (src2 != 0xdeadbeef) + (src3 != 0xdeadbeef);
      U->todo = 1;
      U->INSTTYPE = insn;
      U->prediction_result = (prediction_result == 1);
      if (SafeStride < (1 << 15) - 1)
	SafeStride++;
      if (prediction_result != 2)
	{
	  if (prediction_result)
	    {
	      if (U->predstride)
		if (SafeStride < (1 << 15) - 1)
		  SafeStride += 4 * (1 + (insn == loadInstClass));
	    }
	  else
	    {
	      if (U->predstride)
		SafeStride -= 1024;
	      if (U->predvtage)
		LastMispVT = 0;

	      if (U->predhcvp)
		LastMispHCVP = 0;
               if (U->predvsep)
		LastMispVSEP = 0;


	    }
	}
    }

 if (eligible)
    {
      if (prediction_result == 0)
	NTOTAL++;
      else if (prediction_result < 2)
	NPRED++;
if (U->INSTTYPE== loadInstClass)    {
      if (prediction_result == 0)
	NLOADTOTAL++;
      else if (prediction_result < 2)
	NLOADPRED++;
    }
    }
}

void
beginPredictor (int argc_other, char **argv_other)
{
  (void) argc_other;
  (void) argv_other;

    /** BEGIN FOR HCVP **/
  VP.init (112, 16); // 112, 16
  // A. Seznec seems to be a reasonable trade-off
    /** END FOR HCVP **/
}

void
endPredictor ()
{
 printf ("NPRED: %d ", NPRED);
  printf ("NTOTAL: %d ", NTOTAL);
  printf ("NHCVP: %d ", NHCVP);
  printf ("NVTAGE: %d ", NVTAGE);
  printf ("NSTRIDE: %d ", NSTRIDE);
  printf ("NPREDHCVP: %d ", NPREDHCVP);
  printf ("NPREDVTAGE: %d ", NPREDVTAGE);
  printf ("NPREDSTRIDE: %d ", NPREDSTRIDE);
    printf ("NVSEP: %d ", NVSEP);
    printf ("NPREDVSEP: %d \n", NPREDVSEP);

    printf ("NLOADPRED: %d ", NLOADPRED);
  printf ("NLOADTOTAL: %d ", NLOADTOTAL);
  printf ("NLOADHCVP: %d ", NLOADHCVP);
  printf ("NLOADVTAGE: %d ", NLOADVTAGE);
  printf ("NLOADSTRIDE: %d ", NLOADSTRIDE);
  printf ("NLOADPREDHCVP: %d ", NLOADPREDHCVP);
  printf ("NLOADPREDVTAGE: %d ", NLOADPREDVTAGE);
  printf ("NLOADPREDSTRIDE: %d ", NLOADPREDSTRIDE);
    printf ("NLOADVSEP: %d ", NLOADVSEP);
    printf ("NLOADPREDVSEP: %d \n", NLOADPREDVSEP);
  

}

void
getPredStride (ForUpdate * U, uint64_t & predicted_value, uint64_t seq_no)
{
  bool
    predstride = false;

  uint32_t
    B[NBWAYSTR];
  uint32_t
    TAG[NBWAYSTR];
  uint64_t
    pc = U->pc;
  //use a 3-way skewed-associative structure  for the stride predictor
  for (int i = 0; i < NBWAYSTR; i++)
    {
      //B[i] index in way i ; TAG[i] tag in way i;
      B[i] =
	((((pc) ^ (pc >> (2 * LOGSTR - i)) ^ (pc >> (LOGSTR - i)) ^
	   (pc >> (3 * LOGSTR - i))) * NBWAYSTR) +
	 i) % (NBWAYSTR * (1 << LOGSTR));
      int
	j = (NBWAYSTR - i);
      TAG[i] =
	((pc >> (LOGSTR - j)) ^ (pc >> (2 * LOGSTR - j)) ^
	 (pc >> (3 * LOGSTR - j)) ^ (pc >> (4 * LOGSTR - j))) & ((1 <<
								  TAGWIDTHSTR)
								 - 1);
      U->B[i] = B[i];
      U->TAGSTR[i] = TAG[i];
    }

  int
    STHIT = -1;
  for (int i = 0; i < NBWAYSTR; i++)
    {
      if (STR[B[i]].tag == TAG[i])
	{
	  STHIT = B[i];
	  break;
	}
    }
  U->STHIT = STHIT;
  if (STHIT >= 0)
    if (SafeStride >= 0)
      {				// hit
	uint64_t
	  LastCommitedValue = STR[STHIT].LastValue;

	if (STR[STHIT].conf >= MAXCONFIDSTR / 4)
	  {
	    int
	      inflight = 0;
	    // compute the number of inflight instances of the instruction
                for (uint64_t i = seq_commit + 1; i < seq_no; i++) {
                    inflight += (Update[i & (MAXINFLIGHT - 1)].pc == pc);
                    }
	    predicted_value =
	      (uint64_t) ((int64_t) LastCommitedValue +
			  ((inflight + 1) * ((int64_t) STR[STHIT].Stride)));
	    predstride = true;
          }
        
      }
  

  U->predstride = predstride;
}


// Update of the Stride predictor
// function determining whether to  update or not confidence on a correct prediction
bool
strideupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency,
		  int stride)
{
#define UPDATECONFSTR2 (((!U->prediction_result) || (U->predstride)) && ((random () & ((1 << (NOTLLCMISS + NOTL2MISS + NOTL1MISS + 2*MFASTINST  + 2*(U->INSTTYPE!=loadInstClass))) - 1)) == 0))
#define UPDATECONFSTR1 (abs (stride >= 8) ? (UPDATECONFSTR2 || UPDATECONFSTR2) : (UPDATECONFSTR2))
#define UPDATECONFSTR (abs (stride >= 64) ? (UPDATECONFSTR1 || UPDATECONFSTR1) : (UPDATECONFSTR1))
  return (UPDATECONFSTR &
	  ((abs (stride) > 1) || (U->INSTTYPE != loadInstClass)
	   || ((stride == -1) & ((random () & 1) == 0))
	   || ((stride == 1) & ((random () & 3) == 0))));
//All strides are not equal: the smaller the stride the smaller the benefit (not huge :-))
}

void
UpdateStridePred (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
  int
    B[NBWAYSTR];
  int
    TAG[NBWAYSTR];
  for (int i = 0; i < NBWAYSTR; i++)
    {
      B[i] = U->B[i];
      TAG[i] = U->TAGSTR[i];
    }
  int
    STHIT = -1;

  for (int i = 0; i < NBWAYSTR; i++)
    {
      if (STR[B[i]].tag == TAG[i])
	{
	  STHIT = B[i];
	  break;
	}
    }

  if (STHIT >= 0)
    {
      uint64_t
	LastValue = STR[STHIT].LastValue;
      uint64_t
	Value =
	(uint64_t) ((int64_t) LastValue + (int64_t) STR[STHIT].Stride);
      int64_t
	INTER = abs (2 * ((int64_t) actual_value - (int64_t) LastValue) - 1);

      uint64_t
	stridetoalloc =
	(INTER <
	 (1 << LOGSTRIDE)) ? (uint64_t) ((int64_t) actual_value -
					 (int64_t) LastValue) : 0;

      STR[STHIT].LastValue = actual_value;

      //special case when the stride is not determined
      if (STR[STHIT].NotFirstOcc > 0)
	{
	  if (Value == actual_value)
	    {
	      if (STR[STHIT].conf < MAXCONFIDSTR)
		{
		  if (strideupdateconf
		      (U, actual_value, actual_latency, (int) stridetoalloc))
		    STR[STHIT].conf++;
		}

	      if (STR[STHIT].u < 3)
		if (strideupdateconf
		    (U, actual_value, actual_latency, (int) stridetoalloc))
		  STR[STHIT].u++;
	      if (STR[STHIT].conf >= MAXCONFIDSTR / 4)
		STR[STHIT].u = 3;
	    }
	  else
	    {
		  if (STR[STHIT].conf > (1 << (WIDTHCONFIDSTR - 3)))
		    {
		      STR[STHIT].conf -= (1 << (WIDTHCONFIDSTR - 3));
		    }
		  else
		    {
		      STR[STHIT].conf = 0;
		      STR[STHIT].u = 0;
		    }
STR[STHIT].NotFirstOcc = 0;
	      // this allows to  restart a new sequence with a different   stride

	    }
	}
      else
	{
	  //First occurence
	  //        if (STR[STHIT].NotFirstOcc == 0)
	  if (stridetoalloc != 0)
	    //             if ((stridetoalloc != 0) && (stridetoalloc!=1)  && (((int64_t) stridetoalloc) != -1))
	    // we do not waste the stride predictor storage for stride zero
	    {
	      STR[STHIT].Stride = stridetoalloc;
	    }
	  else
	    {
	      // do not pollute the stride predictor with constant data or with invalid strides
	      STR[STHIT].Stride = 0xffff;
	      STR[STHIT].conf = 0;
	      STR[STHIT].u = 0;
	    }

	  STR[STHIT].NotFirstOcc++;
	}
    }
  else				// the address was not present 
    {
	  int
	    X = random () % NBWAYSTR;
	  bool
	    done = false;
	  // the target entry is not a stride candidate
	  for (int i = 0; i < NBWAYSTR; i++)
	    {
	      STHIT = B[X];
	      if (STR[STHIT].conf == 0)
		{
		  STR[STHIT].conf = 1;	//just to allow not to ejected before testing if possible stride candidate
		  STR[STHIT].u = 0;
		  STR[STHIT].tag = TAG[X];
		  STR[STHIT].Stride = 0;
		  STR[STHIT].NotFirstOcc = 0;
		  STR[STHIT].LastValue = actual_value;
		  done = true;
		  break;
		}
	      X = (X + 1) % NBWAYSTR;
	    }
	  // the target entry has not been useful recently
	  if (!done)
	    for (int i = 0; i < NBWAYSTR; i++)
	      {
		STHIT = B[X];
		if (STR[STHIT].u == 0)
		  {
		    STR[STHIT].conf = 1;
		    STR[STHIT].u = 0;
		    STR[STHIT].tag = TAG[X];
		    STR[STHIT].Stride = 0;
		    STR[STHIT].NotFirstOcc = 0;
		    STR[STHIT].LastValue = actual_value;
		    done = true;
		    break;

		  }
		X = (X + 1) % NBWAYSTR;
	      }
	  //if unable to allocate: age some target entry
	  if (!done)
	    {
	      if ((random () &
		   ((1 <<
		     (2 + 2 * (STR[STHIT].conf > (MAXCONFIDSTR) / 8) +
		      2 * (STR[STHIT].conf >= MAXCONFIDSTR / 4))) - 1)) == 0)
		STR[STHIT].u--;
	    }
    }
}

#ifdef VTAGEON
void
getPredVtage (ForUpdate * U, uint64_t & predicted_value)
{
  bool
    predvtage = false;
  uint64_t
    pc = U->pc;
  uint64_t
    PCindex = ((pc) ^ (pc >> 2) ^ (pc >> 5)) % PREDSIZE;
  uint64_t
    PCbank = (PCindex >> LOGBANK) << LOGBANK;
  for (int i = 1; i <= NHIST; i++)
    {
      U->GI[i] = (gi (i, pc) + (PCbank + (i << LOGBANK))) % PREDSIZE;
      U->GTAG[i] = gtag (i, pc);
    }
  U->GTAG[0] = (pc ^ (pc >> 4) ^ (pc >> TAGWIDTH)) & ((1 << TAGWIDTH) - 1);
  U->GI[0] = PCindex;
  U->HitBank = -1;

  for (int i = NHIST; i >= 0; i--)
    {
      if (Vtage[U->GI[i]].tag == U->GTAG[i])
	{
	  U->HitBank = i;
	  break;
	}
    }


  if (LastMispVT >= DISTMISP)
// when a misprediction is encountered on VTAGE, we do not predict with VTAGE for DISTMISP instructions;
// does not bring significant speed-up, but reduces the misprediction number significantly: mispredictions tend to be clustered       
    if (U->HitBank >= 0)
    {
	predvtage = (Vtage[U->GI[U->HitBank]].conf >= ThresVTAGE);         
	if (predvtage)
	  predicted_value = Vtage[U->GI[U->HitBank]].Val;
        
      }


  U->predvtage = predvtage;
}

bool
vtageupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
#define LOWVAL ((abs (2*((int64_t) actual_value)+1)<(1<<16))+ (actual_value==0))
#define updateconf ((random () & (((1 << (LOWVAL+NOTLLCMISS+2*FASTINST+NOTL2MISS+NOTL1MISS + ((U->INSTTYPE!=loadInstClass) ||NOTL1MISS)       ))- 1)))==0)
  switch (U->INSTTYPE)
    {
    case aluInstClass:
    case fpInstClass:
    case slowAluInstClass:
    case undefInstClass:
    case loadInstClass:
    case storeInstClass:
         //return (updateconf);
         return (updateconf || updateconf || updateconf) ;
      break;
    case uncondIndirectBranchInstClass:
      return (true);
      break;
    default:
      return (false);
    };
}


void
UpdateVtagePred (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
  bool
    ShouldWeAllocate = true;
  if (U->HitBank != -1)
    {
      // there was  an  hitting entry in VTAGE
      uint64_t
	index = U->GI[U->HitBank];
      // the entry may have dissappeared in the interval between the prediction and  the commit

      if (Vtage[index].tag == U->GTAG[U->HitBank])
	{
	  //  update the prediction
	  ShouldWeAllocate = (Vtage[index].Val != actual_value);
	  if (!ShouldWeAllocate)
	    {
	      // the predicted result is correct
                 if (Vtage[index].conf < 32)
		if (vtageupdateconf (U, actual_value, actual_latency))
		  Vtage[index].conf++;
	      if (Vtage[index].u < MAXU)
		Vtage[index].u++;
	    }

	  else
	    {
// misprediction: smart reset
	      if (Vtage[index].conf >= ThresVTAGE)
		{

		  Vtage[index].u = 1;
		  Vtage[index].conf = (MAXCONFID - 4);

		}
	      else
		{
	Vtage[index].Val = actual_value;
        Vtage[index].conf = 0;
	Vtage[index].u = 0;
		}
	    }
	}
    }
    if (ShouldWeAllocate)
      {
	  int
	    DEP = (U->HitBank + 1) + ((random () & 7) == 0);
	  if (U->HitBank == 0)
	    DEP++;
	  if (U->HitBank == -1)
	    {
		DEP = random () & 1;
	    }
	  if (DEP > 1)
	    {
	      for (int i = DEP; i <= NHIST; i++)
		{
		  uint32_t
		    index = U->GI[i];

		  if (Vtage[index].u == 0)
		    {
		      Vtage[index].Val = actual_value;
		      Vtage[index].conf = MAXCONFID- 4;	
		      Vtage[index].tag = U->GTAG[i];
		      break;
		    }
		}
	    }
	  else
	    {
	      for (int j = 0; j <= 1; j++)
		{
		  int
		    i = (j + DEP) & 1;
		  int
		    index = U->GI[i];
		  if (Vtage[index].u == 0)
		    {
		      Vtage[index].Val = actual_value;
		      Vtage[index].conf = MAXCONFID -4;
		      if (U->NbOperand == 0)
			if (U->INSTTYPE == aluInstClass)
			  Vtage[index].conf = MAXCONFID;
		      Vtage[index].tag = U->GTAG[i];
		      break;
		    }
		}
	    }
	}
}
#endif

#ifdef VSEPON
void
getPredVSEP (ForUpdate * U, uint64_t & predicted_value , uint64_t seq_no)
{
  bool
    predvsep = false;
  uint64_t
    pc = U->pc;
  uint64_t     PCindex = ((pc) ^ (pc >> 2) ^ (pc >> 5)) % PREDSIZE;
    uint64_t
    PCbank = (PCindex >> LOGBANK) << LOGBANK;
  for (int i = 1; i <= NHIST; i++)
    {
      U->GI[i] = (gi (i, pc) + (PCbank + (i << LOGBANK))) % PREDSIZE;
      U->GTAG[i] = gtag (i, pc);
    }
  U->GTAG[0] = (pc ^ (pc >> 4) ^ (pc >> TAGWIDTH)) & ((1 << TAGWIDTH) - 1);

    PCindex = ((pc) ^ (pc >> 2) ^ (pc >> 5)) % LCVTSIZE;
  // ugly :-)
  //explodes if more  than 1 M instructions
    while ((LCVT[PCindex].pc!=0) & (LCVT[PCindex].pc!=  U->pc)) PCindex = (PCindex + 67315) % LCVTSIZE;
  LCVT[PCindex].pc=  U->pc;
  U->PCindex = PCindex;

  U->VHitBank = -1;

  for (int i = NHIST; i >= 0; i--)
    {
      if (VSEP[U->GI[i]].tag == U->GTAG[i])
	{
	  U->VHitBank = i;
	  break;
	}
    }


 if (LastMispVSEP >= 128)
// when a misprediction is encountered on VSEP, we do not predict with VSEP for 128 instructions;
// does not bring significant speed-up, but reduces the misprediction number significantly: mispredictions tend to be clustered       
    if (U->VHitBank >= 0)
      {
	predvsep = (VSEP[U->GI[U->VHitBank]].conf>= ThresVSEP);
        U->VSEP_value= LCVT[U->PCindex].Val;
            if (predvsep)
                 predicted_value = U->VSEP_value;
      }
  U->predvsep = predvsep;
}

bool
vsepupdateconf (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
#define LOWVAL ((abs (2*((int64_t) actual_value)+1)<(1<<16))+ (actual_value==0))
#define updateconf ((random () & (((1 << (LOWVAL+NOTLLCMISS+2*FASTINST+NOTL2MISS+NOTL1MISS + ((U->INSTTYPE!=loadInstClass) ||NOTL1MISS)       ))- 1)))==0)
  switch (U->INSTTYPE)
    {
    case aluInstClass:
    case fpInstClass:
    case slowAluInstClass:
    case undefInstClass:
    case loadInstClass:
    case storeInstClass:
//         return (updateconf);
                  return (updateconf||updateconf||updateconf); 
      break;
    case uncondIndirectBranchInstClass:
      return (true);
      break;
    default:
      return (false);
    };
}


void
UpdateVSEPPred (ForUpdate * U, uint64_t actual_value, int actual_latency)
{
  bool
    ShouldWeAllocate = true;
  if (U->VHitBank != -1)
    {
      // there was  an  hitting entry in VSEP
      uint64_t
	index = U->GI[U->VHitBank];
      // the entry may have dissappeared in the interval between the prediction and  the commit
      if (VSEP[index].tag == U->GTAG[U->VHitBank])
	{
	  //  update the prediction
	  ShouldWeAllocate = (U->VSEP_value != actual_value);
	  if (!ShouldWeAllocate)
	    {
	      // the predicted result is correct
	      if (VSEP[index].conf < 32)
		if (vsepupdateconf (U, actual_value, actual_latency))
		  VSEP[index].conf++;
	      if (VSEP[index].u < MAXU)
		VSEP[index].u++;
	    }
	  else
	    {
// misprediction: smart reset
	      if (VSEP[index].conf >= ThresVSEP)
		{
		  VSEP[index].u = 1;
		  VSEP[index].conf = MAXCONFID -(MAXCONFID/4);
		}
	      else
		{
		  VSEP[index].conf = 0;
		  VSEP[index].u = 0;
		}
	    }
	}

    }

LCVT[U->PCindex].Val = actual_value;
LCVT[U->PCindex].pc = U->pc;

    if (ShouldWeAllocate)
      {
	{
	  int
	    DEP = (U->VHitBank + 1) + ((random () & 7) == 0);
	  if (U->VHitBank == 0)
	    DEP++;
	  if (U->VHitBank == -1)
	    {
		DEP = random () & 1;
	    }
	  if (DEP > 1)
	    {
	      for (int i = DEP; i <= NHIST; i++)
		{
		  uint32_t
		    index = U->GI[i];
		  if (VSEP[index].u == 0)
		    {
		      VSEP[index].conf = MAXCONFID - 4;
		      VSEP[index].tag = U->GTAG[i];
		      break;
		    }
		}
	    }

	  else
	    {
	      for (int j = 0; j <= 1; j++)
		{
		  int
		    i = (j + DEP) & 1;

		  int
		    index = U->GI[i];
		  if (VSEP[index].u == 0)
		    {
		      VSEP[index].conf = MAXCONFID - 4;
		      if (U->NbOperand == 0)
			if (U->INSTTYPE == aluInstClass)
			  VSEP[index].conf = MAXCONFID;
		      VSEP[index].tag = U->GTAG[i];
		      break;
		    }
		}
            }
          
	}
      }
}

#endif


