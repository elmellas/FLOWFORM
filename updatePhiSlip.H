Uslip = - tauP
         *
         //(I-tauP*Gr.T())
         //&
        (
          fvc::ddt(U) +
          (U & fvc::grad(U))
        );

phislip = fvc::flux(Uslip);
