      parameter(nlt=33)
      parameter(npr=15)
      data npst,npnd /3,17/
#if RAD2HR
      parameter(nhn=12)
#else
      parameter(nhn=24)
#endif
      parameter(npar=npr)
      common /blam/ lam(nlt)
      common /bnlt/ Fobar(nlt),thray(nlt),oza(nlt),awv(nlt),
     *ao(nlt),aco2(nlt)
      common /bwat/ aw(nlt),bw(nlt)
