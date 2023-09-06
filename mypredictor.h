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

#include <stdint.h>
#include "cvp.h"

#define dfcmh

// ------------------------ fcm -------------------------
#ifdef dfcmh

#include <vector>
#include <deque>
#include <queue>
using namespace std;


struct DFCM_Predictor
{	
	struct FirstLevelEntry
	{
		uint64_t lastValue;
		deque<uint64_t>strideList;
		uint64_t inflight;

		FirstLevelEntry() :lastValue(0), inflight(0), strideList(deque<uint64_t>(4,0)){}
	};

	vector<FirstLevelEntry> Level_1_Table;


	struct SecondLevelEntry
	{
		uint64_t stride;
		uint8_t conf;

		SecondLevelEntry() : stride(0), conf(0) {}
	};
	
	vector<SecondLevelEntry> Level_2_Table;
	
	
	struct InflightInfo
	{
		uint64_t seqNum;
		uint64_t index;
		uint64_t stride_Index;

		InflightInfo() : seqNum(0), index(0) {}
		InflightInfo(uint64_t sn, uint64_t idx, uint64_t str_idx) : seqNum(sn), index(idx), stride_Index(str_idx) {}


	};

	deque<InflightInfo> inflightPreds;
	
	uint64_t globalHistoryRegister;
  	uint64_t pathHistoryRegister;
  	uint64_t addressHistoryRegister;

	uint64_t indexMask;

public:
	
DFCM_Predictor() : globalHistoryRegister(0), pathHistoryRegister(0), addressHistoryRegister(0)
{
	Level_1_Table.resize(4096 * 4);
	Level_2_Table.resize(4096 * 4);
	indexMask = (4096 *4)-1;
	
}
};
	
#endif
// -------------------------------------------------




// Added by A. Seznec
int LastMispHCVP;int LastMispVSEP;
static uint8_t ThresHCVP= 16;
static uint8_t ThresVTAGE= 16;
static uint8_t ThresVSEP= 16;

/*
 *   _    _  _______      _______
 *  | |  | |/ ____\ \    / /  __ \
 *  | |__| | |     \ \  / /| |__) |
 *  |  __  | |      \ \/ / |  ___/
 *  | |  | | |____   \  /  | |
 *  |_|  |_|\_____|   \/   |_|
 *
 * The following code is for HCVP. Above each structure, there is a breakdown of the storage cost.
 * The main code for the predictor is in the ValueCorrelationPredictor class. The total predictor size
 * is given in the comments.
 */
#define DISTMISP 256

#include <vector>
#include <unordered_map>

/*
 * Helper function to generate bit mask (up to 63 bits).
 */
static constexpr uint64_t
makeMask (uint64_t bits)
{
  return ((1ull << bits) - 1);
}

/*
 * Hash function for stride history. Borrowed from DFCM++.
 */
static inline uint64_t
fold (uint64_t result, int bits)
{
  if (bits == 0)
    {
      return result;
    }

  uint64_t val = result;
  val ^= (result >> (64 - bits));
  val &= makeMask (64 - bits);

  return val;
}

/*
 * Structure to maintain meta-data for in-flight instructions.
 *
 * Size:
 *   seq_no: 64 bits
 *   pc: 64 bits
 *   piece: 8 bits
 *   type: 4 bits
 *   eligible: 1 bit
 *   ---
 *   total: 141 bits
 * Does not count toward storage limit.
 */
struct InstInfo
{
  uint64_t seq_no;
  uint64_t pc;
  uint8_t piece;
  InstClass type;
  bool eligible;

    InstInfo (void):InstInfo (0, 0, 0)
  {
  }
  InstInfo (uint64_t seq_no, uint64_t pc, uint64_t piece):seq_no (seq_no),
    pc (pc), piece (piece), type (undefInstClass), eligible (false)
  {
  }
};

/*
 * Maintains a history of branch outcomes as a shift register.
 *
 * Size: 128 bits (only bhrs vector counts toward storage limit)
 */
class BranchHistory
{
private:
  std::vector < uint64_t > bhrs;
  uint64_t size;
  uint64_t bhr_hash;

public:
    BranchHistory (void):bhr_hash (0)
  {
    resize (128);
  }

  void resize (uint64_t size)
  {
    if(size > 1024) return;
    this->size = size;
    bhrs.resize ((size + 63) / 64);
  }

  /*
   * Shift in a single bit outcome for the most recent branch
   */

void shift (uint64_t outcome)
  {
    bhr_hash = 0;
    for (uint64_t i = bhrs.size () - 1; i >= 1; i -= 1)
      {
	bhrs[i] <<= 1;
	bhrs[i] |= (bhrs[i - 1] >> 63);
          }
    bhrs[0] <<= 1;
    bhrs[0] ^= outcome;
    bhrs[bhrs.size () - 1] &= (((1ull <<((this->size-1) % 64)) -1)<<1)+1;

    for (int64_t i = bhrs.size () - 1; i >= 0; i -= 1)bhr_hash ^= bhrs[i];
  
  }

  /*
   * Perform slightly less accurate comparisons on 64-bit hash for simulator performance.
   */
  bool operator== (BranchHistory const &rhs) const
  {
    return hash () == rhs.hash ();
  }
  uint64_t hash (void) const
  {
    return bhr_hash;
  }
};

/*
 * Wraps all the context used to index the predictor tables.
 *
 * Size:
 *   bhr: 128 bits
 *   pc: 64 bits
 *   ---
 *   total: 192 bits
 */
class Context
{
private:
  BranchHistory bhr;
  uint64_t pc;

public:
    Context (void) = default;
    Context (BranchHistory const &bhr, uint64_t pc):bhr (bhr), pc (pc)
  {
  }

  bool valid ()
  {
    return (pc != 0);
  }
  BranchHistory const &getBHR (void) const
  {
    return bhr;
  }
  uint64_t getPC (void) const
  {
    return pc;
  }

  /*
   * Perform slightly less accurate comparisons on 64-bit hash for simulator performance.
   */
  bool operator== (Context const &rhs) const
  {
    return hash () == rhs.hash ();
  }
  uint64_t hash (void) const
  {
    return (pc ^ bhr.hash ());
  }
};

/*
 * Maintains a history of stride values for differential finite context methods.
 *
 * Size: 64 bits
 */
class StrideHistoryEntry
{
private:
  uint64_t stride_hash;

public:
  StrideHistoryEntry (void):stride_hash (0)
  {
  }

  void shift (int64_t stride)
  {
    stride_hash <<= 4;
    stride_hash ^= fold (stride, 16);
  }

  /*
   * Perform slightly less accurate comparisons on 64-bit hash for simulator performance.
   */
  bool operator== (StrideHistoryEntry const &rhs) const
  {
    return hash () == rhs.hash ();
  }
  uint64_t hash (void) const
  {
    return stride_hash;
  }
};

/*
 * An entry that maintains the confidence of each stride that will be added to the base value.
 *
 * Size:
 *   stride: 64 bits
 *   confidence: 4 bits (since counter saturates at 10)
 *   ---
 *   total: 68 bits
 */
class PredictionEntry
{
private:
  int64_t stride;
  uint8_t confidence;

public:
    PredictionEntry (void):PredictionEntry (0, 0)
  {
  }
  PredictionEntry (int64_t stride, uint8_t confidence):stride (stride),
    confidence (confidence)
  {
  }

  int64_t getStride (void) const
  {
    return stride;
  }
  uint8_t getConfidence (void) const
  {
    return confidence;
  }

  void setStride (int64_t stride)
  {
    this->stride = stride;
  }
  void setConfidence (uint8_t confidence)
  {
    this->confidence = confidence;
  }

  /*
   * Perform slightly less accurate comparisons on 64-bit hash for simulator performance.
   */
  bool operator== (PredictionEntry const &rhs) const
  {
    return stride == rhs.stride && confidence == rhs.confidence;
  }
  uint64_t hash (void) const
  {
    return stride ^ confidence;
  }
};

/*
 * Hash functions for use in std::unordered_map.
 */
namespace std
{
  template <> struct hash <StrideHistoryEntry >
  {
    std::size_t operator () (StrideHistoryEntry const &h) const
    {
      return h.hash ();
    }
  };
    template <> struct hash <Context >
  {
    std::size_t operator () (Context const &h) const
    {
      return h.hash ();
    }
  };
    template <> struct hash <BranchHistory >
  {
    std::size_t operator () (BranchHistory const &h) const
    {
      return h.hash ();
    }
  };
    template <> struct hash <PredictionEntry >
  {
    std::size_t operator () (PredictionEntry const &h) const
    {
      return h.hash ();
    }
  };
};

/*
 * Main predictor putting all the pieces together.
 */
class ValueCorrelationPredictor
{
private:
  /*
   * Global BHR.
   *
   * Size: 128 bits
   */
  BranchHistory bhr;

  /*
   * Meta-data that does not count toward storage limit.
   */
  std::unordered_map < Context, int >inflight;
    std::unordered_map < uint64_t, Context > seq_no_to_context;	// take a snapshot of the context used for prediction so that the correct entries are updated

  /*
   * Three primary structures in HCVP.
   *   stride_hist corresponds to Stride History Table (SHT)
   *     Size per entry: 192 bits + 64 bits = 256 bits
   *   prev_val_hist corresponds to Base Value Table (BVT)
   *     Size per entry: 192 bits + 64 bits = 256 bits
   *   predictions corresponds to Value Prediction Table (VPT)
   *     Size per entry: 64 bits + 68 bits = 132 bits
   */
    std::unordered_map < Context, StrideHistoryEntry > stride_hist;
    std::unordered_map < Context, uint64_t > prev_val_hist;
    std::unordered_map < StrideHistoryEntry, PredictionEntry > predictions;

  uint8_t confidence_threshold;

public:
  void init (uint64_t bhr_size, uint8_t confidence_threshold)
  {
    bhr.resize (bhr_size);
    this->confidence_threshold = confidence_threshold;
    ThresHCVP= confidence_threshold;
    
  }

  void shiftBHR (uint64_t outcome)
  {
    bhr.shift (outcome);
  }

  BranchHistory const &getBHR (void) const
  {
    return bhr;
  }

  std::pair < bool, uint64_t > getPrediction (uint64_t seq_no,
					      InstInfo const &info)
  {
    uint64_t pc = info.pc;
    uint8_t piece = info.piece;
    uint64_t pc_idx = (pc << 1) | piece;
    Context context
    {
    bhr, pc_idx};
    seq_no_to_context[seq_no] = context;

    // If the context is currently in flight, do not make a prediction.
    auto inflight_search = inflight.find (context);
    if (inflight_search != inflight.end ())
      {
	inflight_search->second += 1;
	return std::make_pair (false, 0);
      }
    // If the context is not currently in flight, set its counter to 1.
    inflight[context] = 1;

    if (!context.valid ())
      {
	return std::make_pair (false, 0);
      }

    // If the context has not been seen before, do not make a prediction.
    auto stride_hist_search = stride_hist.find (context);
    if (stride_hist_search == stride_hist.end ())
      {
	return std::make_pair (false, 0);
      }

    // If the stride history has not been seen before, do not make a prediction.
    auto prediction_search = predictions.find (stride_hist_search->second);
    if (prediction_search == predictions.end ())
      {
	return std::make_pair (false, 0);
      }

    // Only make a prediction if the confidence is high enough.
//Modified A. Seznec
    if ((prediction_search->second.getConfidence () >= ThresHCVP))

         if (LastMispHCVP>=DISTMISP)
// end
         {
	uint64_t predicted_value =
	  prev_val_hist[context] + prediction_search->second.getStride ();
	return std::make_pair (true, predicted_value);
      }

    return std::make_pair (false, 0);
  }

//    void updatePredictor(uint64_t seq_no, uint64_t actual_value)
  void updatePredictor (uint64_t seq_no, uint64_t actual_value,
			bool UPDATECONFHCVP)
  {
    Context context = seq_no_to_context[seq_no];
    seq_no_to_context.erase (seq_no);

    // This instance of the context is no longer in flight after this point, so update the inflight table.
    auto inflight_search = inflight.find (context);
    if (inflight_search != inflight.end ())
      {
	inflight_search->second -= 1;
	if (inflight_search->second == 0)
	  {
	    inflight.erase (inflight_search);
	  }
      }

    if (!context.valid ())
      {
	return;
      }

    // If this context has not been seen before in the base value table, store the base value.
    auto prev_val_hist_search = prev_val_hist.find (context);
    if (prev_val_hist_search == prev_val_hist.end ())
      {
	prev_val_hist[context] = actual_value;
	return;
      }

    // If the context has been seen before in the base value table, update it.
    int64_t current_delta = actual_value - prev_val_hist_search->second;
    prev_val_hist_search->second = actual_value;

    // If the context has not been seen before in the stride history table, store a new stride history.
    auto stride_hist_search = stride_hist.find (context);
    if (stride_hist_search == stride_hist.end ())
      {
	StrideHistoryEntry stride_hist_entry;
	stride_hist_entry.shift (current_delta);
	stride_hist[context] = stride_hist_entry;
	return;
      }

    // If the context has been seen before in the stride history table, update it.
    StrideHistoryEntry stride = stride_hist_search->second;
    stride_hist_search->second.shift (current_delta);

    // If the stride history has not been seen before, store a new prediction with 0 confidence.
    auto prediction_search = predictions.find (stride);
    if (prediction_search == predictions.end ())
      {
	predictions[stride] = PredictionEntry
	{
//  modified A. Seznec
             current_delta, (uint8_t) (confidence_threshold -4)};
// when a context mispredicts with low confidence it will not be able to reach high confidence unless a big change in the behavior of the application        
	return;
      }

    // If the stride history has been seen before, update the prediction.
    if (prediction_search->second.getStride () ==current_delta)
      {
	// If the stride was the same as last time, increment the confidence.
//A. Seznec
//Confidence increment based on potential benefit
           if (UPDATECONFHCVP)
                if ((abs ((int) current_delta)>=8)  || ((abs ((int) current_delta)>=3) & (random() & 1))|| ((abs ((int) current_delta)<=3) & (((random() & 3)==0) & (current_delta !=0))) || ((random() & 7)==0))
                     prediction_search->second.setConfidence (std::min ( static_cast<uint8_t>  (prediction_search->second.getConfidence () + 1),(uint8_t) 32));//                                                                          confidence_threshold)); 
            
                     
      }
    
    else
	if (prediction_search->second.getConfidence () != 0)
	  {
               prediction_search->second.setConfidence (0);
               
	  }
    else
    {
      // If the stride was different than last time and the confidence was already reset, store the new stride.
      prediction_search->second.setStride (current_delta);
    }
  }
};

/*
 *   ________      ________  _____
 *  |  ____\ \    / /  ____|/ ____|
 *  | |__   \ \  / /| |__  | (___
 *  |  __|   \ \/ / |  __|  \___ \
 *  | |____   \  /  | |____ ____) |
 *  |______|   \/   |______|_____/
 *
 * The following code is a reproduction of the EVES CVP submission with the following modifications:
 *   - Code specific to VTAGE has been removed.
 *   - Parameters for limited storage budget have been removed.
 *   - The ForUpdate structure has had a predhcvp boolean added.
 *   - Formatting has been minorly changed.
 */

/*Copyright (c) <2018>, <2020> INRIA : Institut National de Recherche en Informatique et en Automatique (French National Research Institute for Computer Science and Applied Mathematics)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

/* Same predictor for the 3 tracks, but with different parameters*/

#include <deque>

//#define unlimit
#define k32

#ifdef unlimit
#define UWIDTH 1
#define LOGBANK 20
#define TAGWIDTH 15
#define NBBANK 63

#define NHIST 14
int HL[NHIST + 1] =
  { 0, 0, 1, 3, 7, 15, 31, 47, 63, 95, 127, 191, 255, 383, 511 };
#define LOGSTR 20
#define TAGWIDTHSTR 15
#define LOGSTRIDE 30
#define NBWAYSTR 3
#endif

#ifdef k32

#define UWIDTH 2
#define LOGBANK 9
#define TAGWIDTH 11
#define NBBANK 49

#define NHIST 8
int HL[NHIST + 1] =
  { 0, 0, 3, 7, 15, 31, 63, 90, 127};
#define LOGSTR 4
#define TAGWIDTHSTR 14
#define LOGSTRIDE 20
#define NBWAYSTR 3


#endif

#define WIDTHCONFID 5//4
#define MAXCONFID ((1<< WIDTHCONFID)-1)
#define WIDTHCONFIDSTR 5
#define MAXCONFIDSTR  ((1<< WIDTHCONFIDSTR)-1)
#define MAXU  ((1<< UWIDTH)-1)

#define MINSTRIDE -(1<<(LOGSTRIDE-1))
#define MAXSTRIDE (-MINSTRIDE-1)
#define BANKSIZE (1<<LOGBANK)
#define PREDSIZE (NBBANK*BANKSIZE)

// Global path history

static uint64_t gpath[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

/* using up to 512 bits of path history was found to result in some performance benefit : essentially in the unlimited case. I did not explore longer histories */

static uint64_t gtargeth = 0;
/* history of the targets : limited to  64 bits*/


// The E-Stride predictor
//entry in the stride predictor
struct strdata
{
  uint64_t LastValue;		//64 bits
  uint64_t Stride;		// LOGSTRIDE bits
  uint8_t conf;			// WIDTHCONFIDSTR bits
  uint16_t tag;			//TAGWIDTHSTR bits
  uint16_t NotFirstOcc;		//1 bits
  int u;			// 2 bits
  //67 + LOGSTRIDE + WIDTHCONFIDSTR + TAGWIDTHSTR bits
};
static strdata STR[NBWAYSTR * (1 << LOGSTR)];


static int SafeStride = 0;	// 16 bits
static int LastMispVT = 0;	//8 bits //for tracking the last misprediction on VTAGE

/////////////////////////////////// For E-VTAGE
//VTAGE
struct vtentry
{
  uint64_t Val;		
  uint8_t conf;			//WIDTHCONFID bits
  uint16_t tag;			// TAGWIDTH bits
  uint8_t u;			//2 bits
};

static vtentry Vtage[PREDSIZE];


struct vsepentry
{
  uint8_t conf;			//WIDTHCONFID bits
  uint16_t tag;			// TAGWIDTH bits
  uint8_t u;			//2 bits
};

static vsepentry VSEP[PREDSIZE];

#define LCVTSIZE (1<<20)
struct lcvtentry
{
     uint64_t Val;
       uint64_t pc;
};

static lcvtentry LCVT[LCVTSIZE];

#define  MAXTICK 1024
static int TICK;		//10 bits // for managing replacement on the VTAGE entries



 //index function for VTAGE (use the global path history): just a complex hash function
uint32_t
gi (int i, uint64_t pc)
{
  int hl = (HL[i] < 64) ? (HL[i] % 64) : 64;
  uint64_t inter = (hl < 64) ? (((1 << hl) - 1) & gpath[0]) : gpath[0];
  uint64_t res = 0;
  inter ^= (pc >> (i)) ^ (pc);

  for (int t = 0; t < 8; t++)
    {
      res ^= inter;
      inter ^= ((inter & 15) << 16);
      inter >>= (LOGBANK - ((NHIST - i + LOGBANK - 1) % (LOGBANK - 1)));
    }
  hl = (hl < (HL[NHIST] + 1) / 2) ? hl : ((HL[NHIST] + 1) / 2);

  inter ^= (hl < 64) ? (((1 << hl) - 1) & gtargeth) : gtargeth;
  for (int t = 0; t <= hl / LOGBANK; t++)
    {
      res ^= inter;
      inter ^= ((inter & 15) << 16);
      inter >>= LOGBANK;
    }

  if (HL[i] >= 64)
    {
      int REMAIN = HL[i] - 64;
      hl = REMAIN;
      int PT = 1;

      while (REMAIN > 0)
	{


	  inter ^= ((hl < 64) ? (((1 << hl) - 1) & gpath[PT]) : gpath[PT]);
	  for (int t = 0; t < 8; t++)
	    {
	      res ^= inter;
	      inter ^= ((inter & 15) << 16);

	      inter >>= (LOGBANK -
			 ((NHIST - i + LOGBANK - 1) % (LOGBANK - 1)));

	    }
	  REMAIN = REMAIN - 64;
	  PT++;
	}
    }
  return ((uint32_t) res & (BANKSIZE - 1));
}



//tags for VTAGE: just another complex hash function "orthogonal" to the index function
uint32_t
gtag (int i, uint64_t pc)
{
  int hl = (HL[i] < 64) ? (HL[i] % 64) : 64;
  uint64_t inter = (hl < 64) ? (((1 << hl) - 1) & gpath[0]) : gpath[0];

  uint64_t res = 0;
  inter ^= ((pc >> (i)) ^ (pc >> (5 + i)) ^ (pc));
  for (int t = 0; t < 8; t++)
    {
      res ^= inter;
      inter ^= ((inter & 31) << 14);
      inter >>= (LOGBANK - ((NHIST - i + LOGBANK - 2) % (LOGBANK - 1)));
    }
  hl = (hl < (HL[NHIST] + 1) / 2) ? hl : ((HL[NHIST] + 1) / 2);
  inter ^= ((hl < 64) ? (((1 << hl) - 1) & gtargeth) : gtargeth);
  for (int t = 0; t <= hl / TAGWIDTH; t++)
    {
      res ^= inter;
      inter ^= ((inter & 15) << 16);
      inter >>= TAGWIDTH;
    }

  if (HL[i] >= 64)
    {
      int REMAIN = HL[i] - 64;
      hl = REMAIN;
      int PT = 1;

      while (REMAIN > 0)
	{


	  inter ^= ((hl < 64) ? (((1 << hl) - 1) & gpath[PT]) : gpath[PT]);
	  for (int t = 0; t < 8; t++)
	    {
	      res ^= inter;
	      inter ^= ((inter & 31) << 14);
	      inter >>= (TAGWIDTH - (NHIST - i - 1));


	    }
	  REMAIN = REMAIN - 64;
	  PT++;
	}
    }

  return ((uint32_t) res & ((1 << TAGWIDTH) - 1));
}





#define  MAXTICK 1024


////// for managing speculative state and forwarding information to the back-end
struct ForUpdate
{
  bool predstride;
  bool predhcvp;bool predvsep;
  // Added for HCVP
  bool predvtage;
  bool prediction_result;
     bool predicted;
     
     uint64_t predicted_value;
      uint64_t VSEP_value;
  uint8_t todo;
  uint64_t pc;
  uint32_t GI[NHIST + 1];
  uint32_t GTAG[NHIST + 1];
  uint32_t B[NBWAYSTR];
  uint32_t TAGSTR[NBWAYSTR];
  int STHIT;
  int HitBank;int VHitBank;
  int PCindex;
     
  int8_t INSTTYPE;
  int8_t NbOperand;

};

#define MAXINFLIGHT 256
static ForUpdate Update[MAXINFLIGHT];	// there may be 256 instructions inflight

