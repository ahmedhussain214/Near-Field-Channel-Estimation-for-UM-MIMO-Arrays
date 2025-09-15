This is a code package related to the follow scientific article:

A. Hussain, A. Abdallah and A. M. Eltawil, "[Redefining Polar Boundaries for Near-Field Channel Estimation for Ultra-Massive MIMO
Antenna Array](https://ieeexplore.ieee.org/abstract/document/10988573)," in IEEE Transactions on Wireless Communications, doi: 10.1109/TWC.2025.3564696.

The package contains a simulation environment, based on Matlab, that reproduces all the numerical results and figures in the article. 
We encourage you to also perform reproducible research!

## Abstract
Ultra massive multiple-input-multiple-output (UM-MIMO) technology has emerged as a promising candidate for 6G
networks, offering ultra-high spectral efficiency in wireless systems. The transition to large antenna arrays specially at
high-frequency bands is fundamentally transforming wireless communication from the traditional far-field to the near-field realm. 
This transition poses a distinct challenge in channel estimation due to the associated pilot overhead from large antenna arrays 
and the absence of angular sparsity in near-field spherical wavefronts. 

However, polar-domain sparsity remains achievable, advocating the use of polar codebooks over traditional angle-based ones. 
Nevertheless, the size of the polar codebook presents a significant challenge, necessitating sampling of distance and angle points
across the entire near-field. In this work, we investigate the near-field channel estimation techniques while identifying the 
boundaries of polar domain sparsity with minimal pilot overhead. We propose a novel polar codebook, which leverages our findings 
from sparsity analysis and exploits the beam-focusing properties of the near-field. Unlike existing work, the proposed polar 
codebook design is agnostic to user range information and has considerably reduced dimensions. Capitalizing on this new polar 
codebook, we introduce the beam focused simultaneous orthogonal matching pursuit (BF-SOMP) algorithm for efficient near-field 
channel estimation. To further improve the channel estimation accuracy,we then present a refinement procedure that iterates over 
off-grid angle and range samples to enhance the estimation accuracy.Simulation results demonstrate that the proposed polar
codebook based algorithms outperform contemporary methods in terms of improved normalized mean square error (NMSE) and 
reduced computational complexity. When compared to existing channelestimation methods, the proposed algorithms achieve an
NMSE improvement of 6 − 7 dB at low and high SNR values with 32 pilots, while utilizing a codebook nearly half the size of 
existing ones.

## License and Referencing

If you in any way use this code for research that results in publications, please cite our original article listed above.
