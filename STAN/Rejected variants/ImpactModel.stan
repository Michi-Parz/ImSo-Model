functions {
//////////////////////////////////////////////////////////////////////////////////////////////////////
  // fc
  real fc_fun (real c0, real cL) {
    return sqrt(3) * c0^2 /(pi()*cL * 0.25);
  }
  // f11
  real f11_fun (real c0, real fc, real l1, real l2) {
    return c0^2 * (1/l1^2 + 1/l2^2)/(4*fc);
  }
  // sigma1
  real rad1(real f,real fc) {
    //help2 = if_else(help > 0, 1/sqrt(help), 0.0);
    return 1/sqrt(1-fc/f);
  }
  // sigma2
  real rad2 (real f, real l1, real l2, real c0) {
    return 4*l1*l2*f^2/c0^2;
  }
  // sigma3
  real rad3 (real f, real fc, real l1, real l2, real c0) {
    real u;
    u = 2*(l1+l2);
    return sqrt(pi()*f*u/(16*c0));
  }
  // u/s * c0/f0 * delta1 + delta2
  real rad_delta (real f, real fc, real l1, real l2, real c0) {
    real lambda;
    real delta1;
    real delta2;

    lambda = sqrt(f/fc);
    delta1 = ((1-lambda^2) * log((1+lambda)/(1-lambda)) + 2*lambda)/(4*pi()^2*(1-lambda^2)^1.5);
    delta2 = fc/2 < f ? 0 : (8*c0^2*(1-2*lambda^2))/(fc^2*pi()^4 * l1*l2 * lambda * sqrt(1-lambda^2));

    return (2*(l1+l2)*c0*delta1)/(l1*l2*fc) + delta2;
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////
  // sigma tilde
  real radiation_fac (real f, real fc, real f11, real c0,
                      real l1, real l2, real rad_max, real s3) {
    real s1;
    real s2;
    real sdelta;


    if (f11 <= fc/2){
      if (f>=fc) {
        return fmin(rad1(f, fc),rad_max);
      }


      if (f < fc) {
        sdelta = rad_delta(f, fc, l1,l2,c0);
        s2 = rad2(f, l1, l2, c0);


        if (f < f11 && sdelta > s2) {
            return fmin(s2,rad_max);
        }

        return fmin(sdelta,rad_max);
      }

    }

    //s3 = rad3(f, fc, l1,l2, c0);

    if (f11 > fc/2) {
      if (f < fc){
        s2 = rad2(f, l1, l2, c0);

        if (s2 < s3){
          return fmin(s2,rad_max);
        }
      }
      if (f > fc) {
        s1 = rad1(f, fc);

        if (s1 < s3) {
          return fmin(s1,rad_max);
        }
      }
    }

    return fmin(s3,rad_max);
  }
/////////////////////////////////////////////////////
// Ln
  real impact_sound_L_n (real freq, real mass, real ref_freq,
                              real offset, real fc, real rad_lvl){
    real eta_tot;
    real reverberation_time;
    real Ln;
    
    eta_tot = 0.01 + 0.5/sqrt(freq);
    reverberation_time = 2.2/(freq * eta_tot);
    
    Ln = offset - 30*log10(mass) + 10*log10(reverberation_time) + 10*log10(rad_lvl) + 10*log10(freq/ref_freq);
    
    return(Ln);
  }
  
/////////////////////////////////////////////////////
// Delta L
  real impact_sound_L_delta (real freq, real q, real resonance){
    real L_delta;
    int indicator;
    
    indicator = 0;
    
    if (freq >= resonance) {
      indicator = 1;
    }
    
    
    
    L_delta = q * log10(freq/resonance)*indicator; 
    
    return(L_delta);
  }
  
// Ln': Braucht es nicht als eigene Funktion weil es ja einfach die Differenz aus Ln und Delta L ist...

  
}
