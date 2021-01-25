## tess19

The tess19 software implements the Bayesian nonparametric methods described in Ge et al., Random Tessellation Forests, 2019. This software constructs a
random forest for posterior prediction of categorical data based on real valued predictors. The trees of the random forest are found through SMC inference.
The source, manual and build instructions for tess19 are provided in the Appendices for the paper *Random tessellation forests* (Proceedings of the 33rd Conference on Neural Information Processing Systems), and in this repository.
This software requires the following R packages: optparse, purrr. This software is released under the open source BSD 2-clause license.


## files
- **appendices.pdf**  The Appendices for the paper *Random tessellation forests* . This material includes lemmas and proofs , additional results about the experiments, and also the *tess19* manual.
- **tess19**   The tess19 software.
- **LICENSES**   The  tess19 software license.
 

## Citation
If you use Tess19 in your research, please cite the following publication:

S. Ge, S. Wang, Y.W. Teh, L. Wang, and L.T. Elliott. Random tessellation forests. 2019. Proceedings of the 33rd Conference on Neural Information Processing Systems. Code & Appendices

## LICENSES
tess19 v1.0. Copyright (c) 2019. Shufei Ge and Lloyd T. Elliott.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



## Help
Please send bugs, feature requests and comments to geshf@shanghaitech.edu.cn or shufeig@sfu.ca
