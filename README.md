# Nuclear

The Nuclear R package provides a set of functions to simulate and estimate non-resonant astrophysical S-factors and thermonuclear reaction rates, based on measured nuclear cross sections. The following reactions are currently included: d(p,&gamma;)<sup>3</sup>He, d(d,n)<sup>3</sup>He, d(d,p)<sup>3</sup>H, <sup>3</sup>H(d,n)<sup>4</sup>He, <sup>3</sup>He(<sup>3</sup>He,2p)<sup>4</sup>He, <sup>3</sup>He(&alpha;,&gamma;)<sup>7</sup>Be, <sup>3</sup>He(d,p)<sup>4</sup>He, <span style="color:lightgray"><sup>7</sup>Be(n,p)<sup>7</sup>Li</span>.

#### d(p,&gamma;)<sup>3</sup>He

- [Theoretical S-factor](https://github.com/RafaelSdeSouza/nuclear/blob/master/R/sfactorDpg.R)  
- [Reaction-Rate](https://github.com/RafaelSdeSouza/nuclear/blob/master/R/NumRateDpg.R) 
- [Data](https://github.com/RafaelSdeSouza/nuclear/blob/master/data/dpgdata.RData)   
 

### d(d,n)<sup>3</sup>He

- [Theoretical S-factor](https://github.com/RafaelSdeSouza/nuclear/blob/master/R/sfactorDpg.R)  
- [Reaction-Rate](https://github.com/RafaelSdeSouza/nuclear/blob/master/R/NumRateDpg.R) 
- [Data](https://github.com/RafaelSdeSouza/nuclear/blob/master/data/dpgdata.RData)  

### Installation
devtools::install_github("rafaelsdesouza/nuclear")



### References
[C. Iliadis, K.S. Anderson, A. Coc, F.X. Timmes, and S. Starrfield](http://iopscience.iop.org/article/10.3847/0004-637X/831/1/107/meta), *Bayesian Estimation of Thermonuclear Reaction Rates*, Astrophys. J. 831, 107 (2016).
