# Supervertex clustering for cortical brain parcellation

  Produced and made publicly available by Salim Arslan (contact: name.surname@imperial.ac.uk)
 
  Redistribution and use in source codes, with or without
  modification, are permitted provided that the following conditions
  are met:
 
    - We kindly request that use of this software be cited in publications as:
      S. Arslan, D. Rueckert, "Multi-level parcellation of the cerebral 
      cortex using resting-state fMRI," Proceedings of MICCAI: International 
      Conference on Medical Image Computing and Computer Assisted 
      Intervention. Vol. 9351 of LNCS. Springer, pp. 47-55.
 
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
 
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
 
    - The names of its contributors may be used to endorse or promote 
      products derived from this software without specific prior written 
      permission.
 
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
  IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  

#Notes
Developed in MATLAB R2014a on Ubuntu 14.04

Prerequisites:

(MUST) iso2mesh: a Matlab-based mesh generator 
       http://iso2mesh.sourceforge.net/
       
(COND) Connectome Workbench
       http://www.humanconnectome.org/software/get-connectome-workbench.html

Dataset:

This version has been tested on the Human Connectome Project (HCP) datasets. 
The code has a deep dependency with the HCP data structures. The data can be 
downloaded from https://db.humanconnectome.org. Make sure that you extract 
all the files into the same folder or you may need to modify the main script 
according to your own directory structure.

#Tutorial
Simply run RUN_supervertex_clustering after setting the Path and updating 
the directories from which the data is read.

