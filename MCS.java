import java.io.*;

public class MCS
{    
    static PrintWriter screen = new PrintWriter (System.out, true);

    protected String Element; 
    protected double Z, A, Ro, T, E,Q;
    protected final double mu= 106;
    protected double theta, X0, Xa, beta;
    
    
    public double getThetaL(double inputEnergy, double thickness)
    { 
        double E=inputEnergy; //muon momentum
        double t=thickness;
              
        
        
        double p=Math.sqrt(E*E-mu*mu); //momentum of the muon
        beta=(p/(Math.sqrt((p*p)+(mu*mu))));
        
        X0=716.4*A*(1/(Z*(Z+1)*Math.log((287/Math.sqrt(Z))))); //Average energy loose due to bremsstrahlung and pair production up to a factor of 1/e
        Xa=X0/Ro;
        
        theta=(13.6/(beta*E))*Q*(Math.sqrt(t/Xa))*(1+0.038*Math.log(t/Xa));
        
        return theta;
    }

    public MCS ( String Element, double z, double a, double ro, double t, double e, double q)
    {
        
        Z=z;
        A=a;
        Ro=ro;
        T=t; //thickness of iron
        E=e; //muon momentum
        Q=q;
        //double Q=1; //charge of muon in units of e
        //double theta; // Multiple Coulomb scattering angle

        

        //screen.println(" Name of the element "+Element);

        

        theta=getThetaL(e, t);

        //screen.println(" X0="+X0+"|| Xa="+Xa+"|| Beta="+beta+"|| Theta="+theta);       

        

    }      

}
