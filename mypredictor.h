/*Copyright (c) <2018>, INRIA : Institut National de Recherche en Informatique et en Automatique (French National Research Institute for Computer Science and Applied Mathematics)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

/* Same predictor for the 3 tracks, but with different parameters*/

#include <vector>
#include <deque>

#include <unordered_map>

std::unordered_map<uint64_t, uint8_t> memory;


// Global path history

static uint64_t gpath[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };

/* using up to 512 bits of path history was found to result in some performance benefit : essentially in the unlimited case. I did not explore longer histories */

static uint64_t gtargeth = 0;
/* history of the targets : limited to  64 bits*/


#define  MAXTICK 1024

////// for managing speculative state and forwarding information to the back-end
struct ForUpdate
{
  // we will do prediction or not?
  bool pred;
  // result of prediction
  bool prediction_result;
  // todo = 1 we should do update proccess 
  uint8_t todo;
  uint64_t pc;
  int8_t INSTTYPE;
  int8_t NbOperand;

};

#define MAXINFLIGHT 256
static ForUpdate Update[MAXINFLIGHT];	// there may be 256 instructions inflight
