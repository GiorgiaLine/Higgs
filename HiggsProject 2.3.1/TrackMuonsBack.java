//try to combine all the separate classes into one
//because seriously the number of branches I have is stupid
import java.io.*;
class TrackMuonsBack
    //Class to track the muons back to the origin
{
    //--------Class methods start here--------//
    static PrintWriter screen = new PrintWriter( System.out, true);
    public static double getDelta(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) throws IOException
    {
        //find the vector of the first line, find the vector of the second by plotting from (0,0,0) to (x1,y1,z1)
        //find the angle between two vectors
        
        double Vx_exit = x3-x2;
        double Vy_exit = y3-y2;
        double Vz_exit = z3-z2;
        
        double Vx_orig = x1;
        double Vy_orig = y1;
        double Vz_orig = z1;
        
        double O_dot_E = (Vx_exit*Vx_orig)+(Vy_exit*Vy_orig)+(Vz_exit*Vz_orig);
        
        double modO = Math.sqrt((Vx_orig*Vx_orig)+(Vy_orig*Vy_orig)+(Vz_orig*Vz_orig));
        double modE = Math.sqrt((Vx_exit*Vx_exit)+(Vy_exit*Vy_exit)+(Vz_exit*Vz_exit));
        
        double modPro = modO*modE;
        
        double delta = Math.acos(O_dot_E/modPro);
        
        return delta;
    }
    
    public static double getd(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3) throws IOException
    {
        //find the vector of the first line, find the vector of the second by plotting from (0,0,0) to (x1,y1,z1)
        //find the angle between two vectors
        
        double Vx_exit = x3-x2;
        double Vy_exit = y3-y2;
        
        double Vx_orig = x1;
        double Vy_orig = y1;
        
        double O_dot_E = (Vx_exit*Vx_orig)+(Vy_exit*Vy_orig);
        
        double modO = Math.sqrt((Vx_orig*Vx_orig)+(Vy_orig*Vy_orig));
        double modE = Math.sqrt((Vx_exit*Vx_exit)+(Vy_exit*Vy_exit));
        
        double modPro = modO*modE;
        
        double d = Math.acos(O_dot_E/modPro);
        
        return d;
    }
    
    public static double getBPhi(double x1, double y1, double z1) throws IOException
    {
        //Can I assume the vector we're moving at an angle to is (1,0,0)?
        
        double dot_prod = (x1*1)+(y1*0)+(z1*0);
        double mod_prod = Math.sqrt((x1*x1)+(y1*y1)+(z1*z1));
        
        double BPhi = dot_prod/mod_prod;
        if (BPhi > Math.PI)
        {
            BPhi= BPhi-(Math.PI);
        }
        if (BPhi < -1*(Math.PI))
        {
            BPhi= BPhi+(Math.PI);
        }
        return BPhi;
    }
    public static double getQ(double d) throws IOException
    {
        double mod_d = Math.sqrt(d*d);
        double Q = d/mod_d;
        return Q;
    }
    public static double getLPhi(double BPhi, double d) throws IOException
    {
        double LPhi = BPhi-d;
        if (LPhi > Math.PI)
        {
            LPhi= LPhi-(Math.PI);
        }
        if (LPhi < -1*(Math.PI))
        {
            LPhi= LPhi+(Math.PI);
        }
        return LPhi;
    }
    public static double getP_T(double d, double B, double R) throws IOException
    {
        double mod_d = Math.sqrt(d*d);
        double Q = d/mod_d;
        
        double P_T = (0.3*B*R*Q)/(2*d);
        
        return P_T;
    }
    //alt methods
     public static double getPperp(double x1, double y1)
     {
         double Pperp=Math.sqrt((x1*x1)+(y1*y1));
         return Pperp;
     }
    
     public static double getphi(double x1, double y1)
     {
         double phi = Math.atan(y1/x1);
         if (x1<0 && y1<0)
         phi=phi-(Math.PI);
         if (x1<0 && y1>0)
         phi=phi+(Math.PI);
         return phi;
     }
     public static double getdelta(double Pperp, double Q, double R, double B)
     {
         double delta=(0.3*B*R*Q)/(2*Pperp);
         return delta;
     }
     public static double getPhi(double phi, double delta)
     {
         double Phi=phi+delta;
         return Phi;
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
    
}