"use strict";

// Karl Tarbet   January 2014
// Gradually Varied Flow in Circular Pipe
//

class CircularGVF {
    error_condition: number;
    vel: number;
    E: number;
    Fr: number;// froude number
    fmom: number;
    fr2: number; /* froude number squared */
    y2: number; // alternate momentum depth
    x: number;   /* current x position */
    dx: number;    /* delta x for solution */
    D: number; // diameter of pipe
    L: number;        /* length of canal section */
    Beta: number;
    T: number;   /* top width */
    P: number;   /* perimeter */
    A: number; //area
    c: number;	   /* c in mannings equation  SI c=1.0 ES c=1.486 */
    ddx: number;    /* partial of diameter with respect to x */
    nfirst: number; // used by solver odesol
    dxdid: number;//used by odesol
    constructor(
        public g: number,/* gravity */
        public y: number, 	/* beginning depth for profile */
        public x1: number,   /* beginning x position */
        public x2: number, /* ending x position      */
        public n: number, /* mannings n */
        public So: number, /* Slope of bottom */
        public Q: number,/* flow rate */
        public d1: number,   /* diameter 1 */
        public d2: number,  /* diameter 2 */
        public steps: number,
        public z1: number,  /* elevation */
        public comment: String) {
        this.error_condition = 0;
        this.dx = (this.x2 - this.x1) / (this.steps - 1);
        this.L = Math.abs(this.x2 - this.x1);
        this.ddx = (d2 - d1) / this.L;
        this.nfirst = 1;
        this.dxdid = 0;
    }

    static calc_beta(y: number, D: number) {
        return (Math.acos(1 - 2 * y / D));
    }

    static area(Beta: number, D: number) {
        return (D * D / 4 * (Beta - Math.cos(Beta) * Math.sin(Beta)));
    }

    static perimeter(Beta: number, D: number) {
        return Beta * D;
    }
    static topWidth(Beta: number, D: number) {
        return D * Math.sin(Beta);
    }

    static MomentAhc(B: number, D: number) {
        var hcA = D * D * D / 24 * (3 * Math.sin(B) - 3 * B * Math.cos(B) - Math.pow(Math.sin(B), 3));
        return hcA;
    }
    static momentum(Q: number, g: number, A: number, Ahc: number) {
        return (Ahc + Q * Q / (g * A));

    }


    // determine slope dy/dx at the x position with water depth y
    // slope is called by the solver 'odesol'
    // sets error_condition, and returns slope
    slope(x: number, ygvf: number) {

        //console.log("solving for slope @x=" + x + "  and y= " +ygvf);

        /*******************************/
        /* variables which vary with x */
        /*******************************/
        this.D = this.d1 + this.ddx * Math.abs(x - this.x1);

        //Arguments to acos and acosl must be in
        //the range -1 to 1, or else acos and acosl
        //return NAN and set the global variable errno to
        //EDOM  Domain error

        if (Math.abs(1 - 2 * ygvf / this.D) > 1) {
            this.error_condition = 1;
            return;
        }


        if (ygvf > this.D) {
            this.error_condition = 1;
            return;

        }

        this.Beta = CircularGVF.calc_beta(ygvf, this.D);
        this.A = CircularGVF.area(this.Beta, this.D);
        this.T = CircularGVF.topWidth(this.Beta, this.D);
        this.P = CircularGVF.perimeter(this.Beta, this.D);


        var dax = 0, /* dA/dx    */
            /* height of centroid */
            hc = this.D * this.D * this.D / 24 * (3 * Math.sin(this.Beta) - 3 * this.Beta * Math.cos(this.Beta) - Math.pow(Math.sin(this.Beta), 3)),
            h = this.ddx * 0.00001,
            /*Beta2 */
            Beta2 = CircularGVF.calc_beta(ygvf, this.D + h);

        if (Math.abs(this.ddx) > 0) {
            dax = (CircularGVF.area(Beta2, this.D) - this.A) / 0.00001;
        }

        var Fq = hc / (this.A * this.A) * dax;

        /* slope from mannings eq. */
        var Sf: number = Math.pow(this.Q * this.n / this.c * Math.pow(this.P / this.A, 2 / 3) / this.A, 2);
        this.fr2 = this.Q * this.Q * this.T / (this.g * Math.pow(this.A, 3));

        var dydx: number = (this.So - Sf + this.Q * this.Q * dax / (this.g * Math.pow(this.A, 3)) - Fq) / (1 - this.fr2);

        //console.log("slope = " + dydx);
        return dydx;


    }
    /****************************************************/
    // Returns depth of Critical Flow
    // Fr^2= Q^2 * T /( g * A^3)
    // Solve for Fr^2 =1.0
    /****************************************************/
    static critical(Q: number, D: number, g: number) {

        // Guess 1/2 full
        var Beta: number = 3.14159 / 2.0;

        // Solve for Beta
        var i: number;
        for (i = 0; i < 50; i++) { // max iterations
            var A: number = CircularGVF.area(Beta, D);
            var T: number = D * Math.sin(Beta);
            var Ap: number = CircularGVF.area(Beta + .001, D);
            var F: number = 1. - Q * Q * T / (g * A * A * A);

            var F1: number = 1. - Q * Q * D * Math.sin(Beta + .001) / (g * Ap * Ap * Ap);

            var m: number = (F1 - F) / (.001);
            Beta -= F / m;
            //cout <<" B "<<B;

            if (Beta > 3.14 || Beta < 0) {
                return 0;
                //      break;
                //        nFailed++;
                //       B=nFailed*.05;
                //      cerr <<" "<<i<<" Beta reset ="<<B<<endl;

            }

            if (Math.abs(F) < 0.001) {
                // console.log("F = " + F);
                //  console.log("After " + i + " Iterations");
                break;
            }

        }
        //console.log("After " + i + " Iterations");
        var Y: number = D / 2 * (1 - Math.cos(Beta));
        return (Y);

    }

    /*********************************************************
    Use Mannings Equation to find Normal Depth
    *********************************************************/
    static normal(Q: number, So: number, D: number, n: number, c: number) {
        /* returns -->  y = depth of flow,

        If Normal Depth is not possible
        //??returns --> 0  When Solving for Qn < Q(input) for y=D

        Q = C/n (A**5/3)(S**1/2)/(P**2/3)  in which C=1 for SI units

          and C=1.486 for ES units, n is the channel roughness coefficient,
        A = cross-sectional area, P is the wetted perimeter, and S = slope of
        energy line.  The possible variables and unknowns are:

          D = diameter of circular channels,
          Q = flowrate in cfs or cubic meter per second for ES or SI units,
          n = channel roughness coefficient, and
          S = slope of energy line.
        */

        //Initial guess at Beta
        // From Jeppsons open channel page 47.  Utah State University

        var Qprime: number = n * Q / (c * Math.pow(D, 8 / 3) * Math.sqrt(So));
        var B: number = 2.3286 * Math.pow(Qprime, 0.244);

        //double A, P, F, F1, Y, m;
        var i: number;
        //int nFailed=0;
        //cout << " Beta " << B;
        for (i = 0; i < 50; i++) { // max iterations
            var A: number = CircularGVF.area(B, D);
            var P = B * D;
            //Q = C/n (A**5/3)(S**1/2)/(P**2/3)  in which C=1 for SI units
            var F = Q - c / n * Math.pow(A, 5 / 3) * Math.sqrt(So) / Math.pow(P, 2 / 3);
            var F1 = Q - c / n * Math.pow(CircularGVF.area(B + .001, D), 5 / 3) * Math.sqrt(So) / Math.pow((B + .001) * D, 2 / 3);

            var m = (F1 - F) / (.001);
            B -= F / m;
            //cout << " B " << B;

            if (B > 3.14 || B < 0) {
                return D;
                //      break;
                //        nFailed++;
                //       B=nFailed*.05;
                //      cerr <<" "<<i<<" Beta reset ="<<B<<endl;

            }

            if (Math.abs(F) < 0.001) {
                //console.log("F = " + F);
                // console.log("After " + i + " Iterations");
                break;
            }

        }

        //console.log("After " + i + " Iterations");
        var Y: number = D / 2 * (1 - Math.cos(B));
        return (Y);
    }


    // based on current value of y,
    // find conjugate depth y2
    // use newton method.
    static MomentumDepth2(y: number, D: number, g: number, Q: number) {

        // Trial..
        var B1 = CircularGVF.calc_beta(y, D); // initial Beta
        // initial Momentum
        // console.log("Depth " + y);
        // console.log("Beta case 1 " + B1);
        var M1 = CircularGVF.momentum(Q, g, this.area(B1, D), CircularGVF.MomentAhc(B1, D));
        // console.log( "Target Momentum = "+ M1);

        for (var Bstart = 2.5; Bstart > .1; Bstart -= .02) {
            var B = Bstart;
            for (var counter = 0; counter < 20; counter++) {

                var M = this.momentum(Q, g, this.area(B, D), this.MomentAhc(B, D));
                var Mplus = this.momentum(Q, g, this.area(B + .001, D), this.MomentAhc(B + .001, D));
                var m = (Mplus - M) / (.001);
                B -= (M - M1) / m;
                if (B > 3.14 || B < 0)
                    break;
                if (Math.abs(M - M1) < 0.001)
                    break;
            }

            if (Math.abs(M - M1) < 0.001 && Math.abs(B1 - B) > 0.001) {
                //console.log( "found convergence" );
                //console.log("M = " + M);
                //  console.log( "B = " + B );
                // double tmpy=D/2.*(1-cos(B));
                // double tmpB=Beta(tmpy,D); //
                //  cout <<"tmp B "<<tmpB<<endl;
                return (D / 2. * (1 - Math.cos(B)));
            }
        }
        console.log("failed to find Another B");
        return 0;
    }




    update() {

        var dydx = this.slope(this.x, this.y); // updates A, Beta, etc
        this.vel = this.Q / this.A;
        this.E = this.vel * this.vel / (2 * this.g) + this.y;
        this.Fr = Math.sqrt(this.fr2);
        var Ahc: number = CircularGVF.MomentAhc(this.Beta, this.D);
        this.fmom = CircularGVF.momentum(this.Q, this.g, this.A, Ahc);

    }
 
    /********************************************************/
    // Ordinary differential equations solver Roland Jeppson
    // Returns zero on Error
    // Returns 1 on ok
    //http://web.mit.edu/ehliu/Public/Spring2006/18.304/implementation_bulirsch_stoer.pdf
    //http://books.google.com/books?id=2AWzP4VJ23sC&pg=PA1166&lpg=PA1166&dq=%22roland+jeppson%22+odesol&source=bl&ots=JJoIr9HH2o&sig=oYTo1xdmF0y-ZJqwGtxojWSoJyw&hl=en&sa=X&ei=z5-aUpq1KoX0oASWgIH4DQ&ved=0CCwQ6AEwAA#v=onepage&q&f=true

    // ybeg is input as  y(x1)
    // ybeg is modified at the end of solution to represent y(x2)  and returned
    odesolc(ybeg: number, x1: number, x2: number): number
    // double err, double h1,double hmin,int nstor)
    {
        "use strict";
        // HACK warnings..  the slope function called by this
        // expects a signature slope(number,number).
        // original code was  slope(number,number[]) .. 
        // some code hacked like this ydummy[0] , because  nv=1 (number of values)
        // trying to minimize changes for comparison/debugging to original C code.
        var err = 0.001;
        var h1 = 0.1;
        var hmin = 0.001;
        //     int nstor=1;

        var dxbetw = 0.000001;
        var ngood = 0, nbad = 0, nbetw = 0, ibetw = 0 /*,nstor=1 */, nv = 1;

        var i = 0, j = 0, ji = 0, k = 0, kk = 0, mmn = 0, n = 0, nstp = 0, ntries = 0, ntries1 = 0, nbetw1 = 0;
        var nseq = [2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96];

        var xp: number[] = [0];
        var yp: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var wk1: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var dydx: number[] = [0];
        var xz1: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var fx: number[] = [0, 0, 0, 0, 0, 0, 0];
        var ydumy: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        var dydumy: number[] = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


        var bbb, binter, cc1c, dumyx = 0, dvv, dxsta, errorm, h, h2, hdid, he;

        var x = 0, x1x = 0, xbetw = 0, xest = 0, xsign = 0, yente = 0, ystor = 0;


        ntries = 7;
        ntries1 = ntries - 1;
        nbetw1 = nbetw - 1;
        xsign = x2 - x1;
        x = x1;
        if (this.nfirst === 1) {
            if (xsign < 0) dxsta = -Math.abs(h1);
            else dxsta = Math.abs(h1);

        }
        else {
            if (xsign < 0) dxsta = -Math.abs(this.dxdid);
            else dxsta = Math.abs(this.dxdid);

        }

        ngood = 0;
        nbad = 0;
        ibetw = 0;

        for (i = 0; i < nv; i++) {
            wk1[i] = ybeg; //ybeg[i];
        }
        xbetw = x - 2 * dxbetw;

        for (nstp = 1; nstp < 10000; nstp++) {

            dydx[0] = this.slope(x, wk1[0]);

            if (this.error_condition === 1) {
                return 0;
            }


            for (i = 0; i < nv; i++)
                wk1[nv + i] = Math.abs(wk1[i]) + Math.abs(dxsta * dydx[i]) + 1.0e-25;
            if (!(nbetw < 1)) {
                if (!(Math.abs(x - xbetw) <= Math.abs(dxbetw))) {
                    if (!(ibetw >= nbetw1)) {
                        xp[ibetw] = x;
                        for (i = 0; i < nv; i++)
                            yp[nv * ibetw + i] = wk1[i];
                        ibetw++;
                        xbetw = x;
                    }
                }
            }

            if ((x + dxsta - x2) * (x + dxsta - x1) > 0.) dxsta = x2 - x;
            h = dxsta;
            xbetw = x;

            for (i = 0; i < nv; i++) {
                wk1[3 * nv + i] = wk1[i];
                wk1[4 * nv + i] = dydx[i];
            }
            /*L60: */
            do {
                for (i = 0; i < 11; i++) {
                    he = h / nseq[i];

                    for (ji = 0; ji < nv; ji++) {
                        wk1[5 * nv + ji] = wk1[3 * nv + ji];
                        ydumy[ji] = wk1[3 * nv + ji] + wk1[4 * nv + ji] * he;
                    }
                    x1x = xbetw + he;
                    dydumy[0] = this.slope(x1x, ydumy[0]);
                    if (this.error_condition === 1) {
                        return 0;
                    }


                    h2 = 2. * he;
                    mmn = Math.round(nseq[i]);
                    for (n = 1; n < mmn; n++) {

                        for (ji = 0; ji < nv; ji++) {
                            ystor = wk1[5 * nv + ji] + h2 * dydumy[ji];
                            wk1[5 * nv + ji] = ydumy[ji];
                            ydumy[ji] = ystor;
                        }
                        x1x = x1x + he;
                        dydumy[0] = this.slope(x1x, ydumy[0]);
                        if (this.error_condition === 1) {
                            return 0;
                        }



                    }
                    for (ji = 0; ji < nv; ji++)
                        dydumy[ji] = (wk1[5 * nv + ji] + ydumy[ji] + he * dydumy[ji]) / 2;
                    xest = h / nseq[i];
                    xest = xest * xest;
                    xz1[i] = xest;
                    if (i == 0) {
                        for (j = 0; j < nv; j++) {
                            wk1[j] = dydumy[j];
                            wk1[6 * nv + j] = dydumy[j];
                            wk1[2 * nv + j] = dydumy[j];
                        }
                    }
                    else {
                        if (i < ntries1) mmn = i + 1;
                        else mmn = ntries;
                        for (k = 1; k < mmn; k++) fx[k] = xz1[i - k] / xest;
                        for (j = 0; j < nv; j++) {
                            yente = dydumy[j];
                            dvv = wk1[6 * nv + j];
                            cc1c = yente;
                            wk1[6 * nv + j] = yente;
                            for (k = 1; k < mmn; k++) {
                                kk = k + 6;
                                binter = fx[k] * dvv;
                                bbb = binter - cc1c;
                                if (bbb != 0.) {
                                    bbb = (cc1c - dvv) / bbb;
                                    dumyx = cc1c * bbb;
                                    cc1c = binter * bbb;
                                }
                                else dumyx = dvv;
                                dvv = wk1[nv * kk + j];
                                wk1[nv * kk + j] = dumyx;
                                yente = yente + dumyx;
                            }
                            wk1[2 * nv + j] = dumyx;
                            wk1[j] = yente;
                        }
                    }
                    errorm = 0;
                    for (j = 0; j < nv; j++)
                        if (Math.abs(wk1[2 * nv + j] / (wk1[nv + j])) > errorm)
                            errorm = Math.abs(wk1[2 * nv + j] / (wk1[nv + j]));
                    if (errorm < err) {
                        x = x + h;
                        hdid = h;
                        if (i == ntries1) this.dxdid = .95 * h;
                        else if (i == ntries1 - 1) this.dxdid = 1.2 * h;
                        else {
                            this.dxdid = h * nseq[ntries1 - 1] / nseq[i];
                        }
                        break;
                    }
                } /*for i*/

                if (!(errorm < err)) {
                    h = 0.0625 * h;
                    if (h < 1.0e-8) {
                        console.log("\nincr.=0");
                        this.error_condition = 1;
                        return (0);
                    }
                }

            } while (errorm > err);
            if (hdid == dxsta) ngood++;
            else nbad++;
            if (!((x - x2) * xsign < 0.)) {
                for (i = 0; i < nv; i++)
                    ybeg = wk1[i];
                if (nbetw <= 0) {
                    this.nfirst = 0;
                    return ybeg;
                }
                xp[ibetw] = x;
                for (i = 0; i < nv; i++)
                    yp[nv * ibetw + i] = wk1[i];
                ibetw++;
                this.nfirst = 0;
                return ybeg;
            }

            if (Math.abs(this.dxdid) < hmin) {
                console.log("\nincrement<min.");
                this.error_condition = 1;
                return (0);

            }

            dxsta = this.dxdid;

        } /*for nstp*/

        console.log("\nallowable increments exceeded");
        this.error_condition = 1;
        return (0);

    } /* odesolc */


    printStartingCondition() {
        console.log("#gradually varied open channel flow");
        console.log("#y=" + this.y);
        console.log("#x1= " + this.x1);
        console.log("#x2= " + this.x2);
        console.log("#dx= " + this.dx);
        console.log("#steps=" + this.steps);
        console.log("#n= " + this.n);
        console.log("#So= " + this.So);
        console.log("#Q= " + this.Q);
        console.log("#z1= " + this.z1);
        console.log("#d1= " + this.d1);
        console.log("#d2= " + this.d2);
        console.log("#ddx= " + this.ddx);
        console.log("#g  = " + this.g);
        console.log("\n X       Y      y/D      yn   Beta  Area    E      M       V       Fr    Z     H        D      y2     abs(y2-y)  yc");
        console.log(this.comment);
    }

    updateResults(results,ygvf,yn,y2,z,yc) {

            console.log(this.x + " " + ygvf + " " + (ygvf / this.D) + " " + yn + " " + this.Beta + " " + this.A + " " +
                this.E + " " + this.fmom + " " + this.vel + " " + this.Fr + " " + z + " " + (z + this.E) +
                " " + this.D + " " + y2 + " " + Math.abs(y2 - this.y) + " " + yc);
            var row = [this.x, ygvf, (ygvf / this.D), yn, this.Beta, this.A,
                this.E, this.fmom, this.vel, this.Fr, z, (z + this.E),
                this.D, y2, Math.abs(y2 - this.y), yc];

            results.resultTable.push(row);
        }

    solve() {
        var ygvf: number;
        var xz: number;	   /* new positon of x   */
        var z: number;
        var rval = [];
        this.x = this.x1;
        this.c = 1; //s.i. units
        if (this.g > 30)
            this.c = 1.486; // english units
        this.D = this.d1;
        ygvf = this.y;
      
        this.printStartingCondition();
        

        var results = { resultTable: [], columnHeader: [] };

        results.columnHeader = ["X", "Y", "y/D", "yn", "Beta", "Area", "E", "M", "V", "Fr", "Z", "H", "D", "y2", "abs(y2-y)", "yc"];
        
        z = this.z1;
       
        /* GVF profile calculations between  */
        for (var istep = 0; istep < this.steps; istep++) {

           var y2 = CircularGVF.MomentumDepth2(this.y, this.D, this.g, this.Q);
           var yn = CircularGVF.normal(this.Q, this.So, this.D, this.n, this.c);
            var yc = CircularGVF.critical(this.Q, this.D, this.g);
            if (this.error_condition === 1) {
                return 0;
            }
            this.update(); // update some variables  
            this.updateResults(results, ygvf, yn, y2, z,yc);

            if (istep === this.steps) {
                break;
            }

            xz = this.x + this.dx;
            console.log("xz=" + xz);
            z = (this.x1 - xz) * this.So + this.z1;
            if (this.error_condition === 0) {
                ygvf = this.odesolc(ygvf, this.x, xz);
                console.log("ygvf= " + ygvf);
            }
            else {
                return 0;
            }

            this.x = xz;
            this.y = ygvf;
            this.update();
           
            if ((ygvf / yn) < 0.1)
                this.error_condition = 1;

            //if( fabs( Fr -1.0) < 0.03)
            //error_condition=1;


        }// while (rep);
        return results;
    }
  
}



//window.onload = () => {
//    var el = document.getElementById('content');
//    var greeter = new Greeter(el);

//    greeter.start();
//};
//
//console.log("hi");
//    colHeaders: ["DIABEGIN", "DIAEND", "SLOPE", "N", "YBEGIN", "XBEGIN", "XEND", "STEPS", "COMMENT"],
// ["7.5", "7.5", "0.1", "0.024", "3.6", "0", "20", "10", "sect 1"],
//var g = 32.2;
//var diameter1 = 7.5;
//var diameter2 = 7.5;
//var slope = .1;
//var n = .024;
//var z = 660.1;
//var y1 = 3.6;
//var gvf = new CircularGVF(g, y1, 0, 20, n, slope, 305.6, diameter1, diameter2, 10, z, "sec1");
//gvf.solve();