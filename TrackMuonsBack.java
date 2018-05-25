//try to combine all the separate classes into one
//because seriously the number of branches I have is stupid
import java.io.*;
class TrackMuonsBack
//Class to track the muons back to the origin
{
    //--------Class methods start here--------//
    static PrintWriter screen = new PrintWriter( System.out, true);

    public static double getQ(double d) throws IOException
    {
        double Q = 1;
        if (d<0)
        {
            Q=-1;
        }
        return Q;
    }
    //alt methods
    public static double getPperp(double B, double R, double Q, double d)
    {
        double Pperp=(0.3*B*R*Q)/(2.0*d);
        return Pperp;
    }

    public static double getPhi(double grad, double x0, double y0)
    {
        double Phi = Math.atan(grad);
        /*if (Phi > 2*Math.PI)
        {
            Phi= Phi-(2*Math.PI);
        }
        if (Phi < 0)
        {
            Phi= Phi+(2*Math.PI);
        }*/
        if(x0 < 0.0 && y0 > 0.0) 
        {
            Phi = Phi + Math.PI;
        }
        if(x0 < 0.0 && y0 < 0.0) 
        {
            Phi = Phi + Math.PI;
        }
        if(x0 > 0.0 && y0 < 0.0) 
        {
            Phi = Phi + (2.0*Math.PI);
        }

        return Phi;
    }

    public static double getLPhi(double BPhi, double delta)
    {
        double LPhi = BPhi - delta;
        if (LPhi > 2*Math.PI)
        {
            LPhi= LPhi-(2*Math.PI);
        }
        if (LPhi < 0)
        {
            LPhi= LPhi+(2*Math.PI);
        }
        return LPhi;
    }

    public static double getdelta(double Phi, double BPhi)
    {
        double delta=Phi-BPhi;
        return delta;
    }

    public static double getBPhi(double x0, double y0)
    {
        double BPhi=Math.atan2(y0, x0);
        if (BPhi > 2*Math.PI)
        {
            BPhi= BPhi-(2*Math.PI);
        }
        if (BPhi < 0)
        {
            BPhi= BPhi+(2*Math.PI);
        }
        /*if(x0 < 0.0 && y0 > 0.0) 
        {
            BPhi = BPhi + Math.PI;
        }
        if(x0 < 0.0 && y0 < 0.0) 
        {
            BPhi = BPhi + Math.PI;
        }
        if(x0 > 0.0 && y0 < 0.0) 
        {
            BPhi = BPhi + (2.0*Math.PI);
        }*/
        return BPhi;
    }

    //Get origional particle specs

    public static double getOriginalMass(double P_T1, double P_T2, double P_T3, double P_T4,
    double Pz1, double Pz2, double Pz3, double Pz4,
    double E1, double E2, double E3, double E4) throws IOException
    //Get Es and Pzs form dataset, P_Ts from get P_T method
    {
        double ETot = E1+E2+E3+E4;
        double PzTot = Pz1+Pz2+Pz3+Pz4;
        double P_TTot = P_T1+P_T2+P_T3+P_T4;
        double Abs_P = Math.sqrt((P_TTot*P_TTot)+(PzTot*PzTot));

        double M_0 = Math.sqrt((ETot*ETot)-(Abs_P*Abs_P));

        return M_0;
    }

    public static double getInvMass( double Prec[][]) throws IOException
    //Get Es and Pzs form dataset, P_Ts from get P_T method
    {
        double Sumpx = 0.0;
        double Sumpy = 0.0;
        double Sumpz = 0.0;
        double SumE  = 0.0;
        for(int i=1; i<=4; i++) {
            Sumpx = Sumpx + Prec[0][i-1];
            Sumpy = Sumpy + Prec[1][i-1];
            Sumpz = Sumpz + Prec[2][i-1];
            SumE  = SumE  + Prec[5][i-1]; 
        }
        double Sumpx2 = Sumpx*Sumpx;
        double Sumpy2 = Sumpy*Sumpy;
        double Sumpz2 = Sumpz*Sumpz;
        double SumE2  = SumE*SumE;
        double InvMassq = SumE2 - Sumpx2 - Sumpy2 - Sumpz2;
        double InvMass = 0.0;
        if(InvMassq > 0.0) {
            InvMass = Math.sqrt(InvMassq);
        }

        return InvMass;
    }

    public static double getAbsMomentum(double P_T1, double P_T2, double P_T3, double P_T4,
    double Pz1, double Pz2, double Pz3, double Pz4) throws IOException
    //Get P's from main method. Pz from dataset, P_T from getP_T method.
    {
        double PzTot = Pz1+Pz2+Pz3+Pz4;
        double P_TTot = P_T1+P_T2+P_T3+P_T4;

        double abs_P = Math.sqrt((P_TTot*P_TTot)+(PzTot*PzTot));

        return abs_P;
    }

    public static double getOriginalQ(double q1, double q2, double q3, double q4) throws IOException
    //Get Qs for each muon im main method, getQ.set1, get\Q.set2,... ect
    {
        double charge = q1+q2+q3+q4;

        return charge;
    }

    //Line fit

    public static double[] lineFit(double x[], double y[], int n) throws IOException
    //line fit code
    {
        double [] Ans;
        Ans = new double [2];

        if(n < 3) 
        {
            Ans[0] = 0.0;
            Ans[1] = 0.0;
        }
        double Count = 0.0;
        double Sumx  = 0.0;
        double Sumy  = 0.0;
        double Sumxy = 0.0;
        double Sumxx = 0.0;
        double Sumyy = 0.0;
        for(int j=1; j<=(n-1); j++) 
        {
            if(y[j-1] != 0.0) 
            {
                Sumx = Sumx + x[j-1];
                Sumy = Sumy + y[j-1];
                Count= Count+ 1.0;
            }
        }
        if(Count <= 1.0) {
            Ans[0] = 0.0;
            Ans[1] = 0.0;
        }
        double Ymed = Sumy/Count;
        double Xmed = Sumx/Count;
        for(int j=1; j<(n-1); j++) {
            if(y[j-1] != 0.0) {
                double Scartx = x[j-1] - Xmed;
                double Scarty = y[j-1] - Ymed;
                Sumxy  = Sumxy + Scartx*Scarty;
                Sumxx  = Sumxx + Scartx*Scartx;
                Sumyy  = Sumyy + Scarty*Scarty;
            }
        }
        // -----------------------------------------------
        // Fit Parameters:
        // -----------------------------------------------
        if(Sumxx == 0.0) {
            Ans[0] = 0.0;
            Ans[1] = 0.0;
        }
        double A = Sumxy/Sumxx;
        double B = Ymed - A*Xmed;
        double E = 0.0;
        if(Count >= 3.0) {
            E = (Sumyy - Sumxy*A)/(Count-2.0);
        }
        Ans[0] = A;
        Ans[1] = B;
        return Ans;
    }
    // --------------------------------------------------------------------
    //float xxx[3];
    //float yyy[3];
    //float Ans[2];
    // --------------------------------------------------------------------
    // Test Linefit:
    // --------------------------------------------------------------------
    // xxx[0] = 1.0;
    // yyy[0] = 4.0*xxx[0] + 3.5;
    // xxx[1] = 2.0;
    // yyy[1] = 4.0*xxx[1] + 3.5;
    // xxx[2] = 3.0;
    // yyy[2] = 4.0*xxx[2] + 3.5;
    // int n3 = 3;
    // Linefit(xxx,yyy,n3,Ans);
    // cout<<" Gradient "<<Ans[0]<<endl;
    // cout<<" Intercept "<<Ans[1]<<endl;

}
