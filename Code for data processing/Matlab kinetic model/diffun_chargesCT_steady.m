function    dcdt = diffun(t,c,k,i0)
            S=c(1); CT=c(2); C=c(3);
            kradnr =k(1); 
            kct    =k(2);  
            kcr    =k(3);
            knrct  =k(4);
            kb     =k(5);
            kcs    =k(6);
            kenc   =k(7);
            knrc   =k(8);
            dcdt(1,1) = i0-(kradnr+kb*(S+CT)+kct)*S+kcr*CT;
            dcdt(2,1) = kct*S-(kcr+kcs)*CT+kenc*C^2;
            dcdt(3,1) = kcs*CT-kenc*C^2-knrc*C;
end
