import java.io.*;
import java.util.*;
public class Import3
{
    static PrintWriter screen = new PrintWriter( System.out, true);

    static String [] club = new String[20]; 
    static int Events = 2022; 
    static int entries_per_particle = 5;
    static int muons_per_event = 4;
    static double [] [] [] Finaldet1;
    static double [] [] [] Finaldet2;
    static double [] [] [] Finaldet3;

    static Histogram HiggsMass = new Histogram("HiggsMass", 1000, 0, 1000);

    public static void main () throws IOException
    {
        BufferedReader dataBR = new BufferedReader(new FileReader(new File("C:\\Users\\Giorgia\\Documents\\Data.csv")));
        String line = "";
        double E =0;    
        double Px = 0;  
        double Py = 0; 
        double Pz = 0;
        double Q = 0;   
        int total = 0;

        double B = 4; //B Field
        double R = 1.2;// Radius of the cavity

        Finaldet1 = new double [2][3][getNumber()];
        Finaldet2 = new double [2][3][getNumber()];
        Finaldet3 = new double [2][3][getNumber()];

        ArrayList<String[]> dataArr = new ArrayList<String[]>();
        ArrayList<double[][]> MuonData = new ArrayList<double[][]>();

        double [] [] MuonStore  = new double [Events*muons_per_event+1] [entries_per_particle];  
        while ((line = dataBR.readLine()) != null) 
        { // Read a single line from the file until there are no more lines to read
            String [] value = line.split("\\,");   

            double [][] muon_data_line = new double[muons_per_event][entries_per_particle];

            for (int i = 0; i < muons_per_event; i++)
            {
                for (int j = 0; j < entries_per_particle; j++)
                {
                    muon_data_line[i][j] = Double.valueOf(value[entries_per_particle*i + j]);
                }
            }

            MuonData.add(muon_data_line);

        }

        for (int k = 0; k < MuonData.size(); k++) 
        {
            for (int i = 0; i < muons_per_event; i++)
            {
                for (int j = 0; j < entries_per_particle; j++)
                {
                    //System.out.println((MuonData.get(k))[i][j]);
                }
            }
        }

        for (int k = 0; k < MuonData.size(); k++) 
        {
            for (int i = 0; i < muons_per_event; i++)
            {
                total = total+1;
                //E  = (MuonData.get(k))[i][0];
                //Px = (MuonData.get(k))[i][1];
                //Py = (MuonData.get(k))[i][2];
                //Pz = (MuonData.get(k))[i][3];
                //Q  = (MuonData.get(k))[i][4];

                //double [] [] [] E  = new double [sizek][sizei][size0];  
                MuonStore[total][0] = MuonData.get(k) [i] [0];  //store energy
                MuonStore[total][1] = MuonData.get(k) [i] [1];  //store px
                MuonStore[total][2] = MuonData.get(k) [i] [2];  //store py
                MuonStore[total][3] = MuonData.get(k) [i] [3];  //store pz
                MuonStore[total][4] = MuonData.get(k) [i] [4];  //store Q

                double p = Math.sqrt(Px*Px + Py*Py + Pz*Pz);
                //System.out.println(" Event number: " + k + " muon number: " + i + " Px " + Px);
            }
        }
        //System.out.println(" Event number Num  5, Muon 0, energy" + MuonStore[20][0]);
        //System.out.println(" Event number Num  6, Muon 1, Px" + MuonStore[25][1]);

        //TrackMuon.TrackerInfo();
        TrackMuon.Tracking(MuonStore); 
        Finaldet1 = TrackMuon.GetFinald1();
        int muonsPerEvent = 4;
        double [][]Prec;
        Prec = new double [6][muonsPerEvent];
        double []InvMass_array;
        InvMass_array = new double [MuonData.size()];
        //double [][][][]Finaldet;
        //Finaldet = new double [2][3][MuonData.size()][muonsPerEvent];

        /*for (int j = 0; j < getNumber()-1; j++) 
        {
        int k = j;
        int i=0;
        for (i=0; (k-4)> 0;) 
        {
        i=i+1;
        k = j-(i*4);
        }
        Finaldet [0][0][i][j-(i*4)]=Finaldet1[0][0][j];
        Finaldet [1][0][i][j-(i*4)]=Finaldet1[1][0][j];
        Finaldet [0][1][i][j-(i*4)]=Finaldet1[0][1][j];
        Finaldet [1][1][i][j-(i*4)]=Finaldet1[1][1][j];
        Finaldet [0][2][i][j-(i*4)]=Finaldet1[0][2][j];
        Finaldet [1][2][i][j-(i*4)]=Finaldet1[1][2][j];

        }*/

        for (int i = 0; i < MuonData.size(); i++) 
        {
            for (int k = 0; k < muonsPerEvent; k++) 
            { 
                //coordinates of detector hits
                double x1= (Finaldet1 [0][0][(i*4)+k])/100;
                double y1= (Finaldet1 [1][0][(i*4)+k])/100;
                //double z1= MuonStore [i][3];

                double x2= (Finaldet1 [0][1][(i*4)+k])/100;
                double y2= (Finaldet1 [1][1][(i*4)+k])/100;
                //double z2= MuonStore [i][3];

                double x3= (Finaldet1 [0][2][(i*4)+k])/100;
                double y3= (Finaldet1 [1][2][(i*4)+k])/100;
                //double z3= MuonStore [i][3];

                double []x_array;
                x_array = new double [3];
                x_array [0] = x1;
                x_array [1] = x2;
                x_array [2] = x3;

                double []y_array;
                y_array = new double [3];
                y_array [0] = y1;
                y_array [1] = y2;
                y_array [2] = y3;

                double [] Ans;
                Ans = new double [2];
                Ans = TrackMuonsBack.lineFit( x_array, y_array, muonsPerEvent);

                double Grad = Ans [0];
                double Inter= Ans [1];

                double xp1 = (-2.0*Grad*Inter + Math.sqrt(4.0*Grad*Grad*Inter*Inter-4.0*(1.0+Grad*Grad)*(Inter*Inter-R*R)))/(2.0*(1.0+Grad*Grad));
                double xp2 = (-2.0*Grad*Inter - Math.sqrt(4.0*Grad*Grad*Inter*Inter-4.0*(1.0+Grad*Grad)*(Inter*Inter-R*R)))/(2.0*(1.0+Grad*Grad));

                double yp1 = Grad*xp1 + Inter;
                double yp2 = Grad*xp2 + Inter;

                double Dist1 = Math.sqrt((x_array[0]-xp1)*(x_array[0]-xp1)+(y_array[0]-yp1)*(y_array[0]-yp1));
                double Dist2 = Math.sqrt((x_array[0]-xp2)*(x_array[0]-xp2)+(y_array[0]-yp2)*(y_array[0]-yp2));
                double x0  = xp1;
                double y0  = yp1;
                if(Dist2 < Dist1) 
                {
                    x0 = xp2;
                    y0 = yp2;
                }

                //double q= MuonStore[i][4];//charge of the muon

                double Phi = TrackMuonsBack.getPhi(Grad, x0, y0);
                double BPhi = TrackMuonsBack.getBPhi(x0, y0);
                double d = TrackMuonsBack.getdelta(Phi, BPhi);
                double LPhi = TrackMuonsBack.getLPhi(BPhi, d);
                double q = TrackMuonsBack.getQ(d);
                double Pperp = TrackMuonsBack.getPperp(B, R, q, d);

                Prec [0][k] = Pperp*Math.cos(LPhi);
                Prec [1][k] = Pperp*Math.sin(LPhi);
                //Prec [2][k] = MuonStore[k][3];
                Prec [2][k] = MuonData.get(i) [k] [3];
                Prec [3][k] = Math.sqrt(Prec[0][k]*Prec[0][k]+Prec[1][k]*Prec[1][k]+Prec[2][k]*Prec[2][k]);
                Prec [4][k] = q;
                Prec [5][k] = Math.sqrt(Prec[3][k]*Prec[3][k] + 0.1057*0.1057);

                /*double [][]PT_array;
                PT_array = new double [ MuonData.size()][muonsPerEvent];

                double [][]Pz_array;
                Pz_array = new double [ MuonData.size()][muonsPerEvent];

                double [][]E_array;
                E_array = new double [ MuonData.size()][muonsPerEvent];

                for (int j = 0; j < muonsPerEvent; j++) 
                {
                PT_array[i][j] = TrackMuonsBack.getPperp(B, R, q, d);

                Pz_array[i][j] = MuonStore[i][j];

                E_array[i][j] = MuonStore[i][j];

                //screen.println("Pz= " +Pz_array[i][j]+ ", PT = " +PT_array[i][j]);
                //screen.println("E= " +E_array[i][j]);

                }/*
                /*double [][]Q_array;
                Q_array = new double [2022][4];
                Q_array[i][0] = TrackMuonsBack.getQ(d);
                Q_array[i][1] = TrackMuonsBack.getQ(d);
                Q_array[i][2] = TrackMuonsBack.getQ(d);
                Q_array[i][3] = TrackMuonsBack.getQ(d);*/

                //double parentMass = TrackMuonsBack.getOriginalMass(PT_array[i][0],PT_array[i][1],PT_array[i][2],PT_array[i][3],
                //Pz_array[i][0],Pz_array[i][1],Pz_array[i][2],Pz_array[i][3],
                //E_array[i][0],E_array[i][1],E_array[i][2],E_array[i][3]);

                //Q = TrackMuonsBack.getQ(d);
                //screen.println("Gradient = " +Grad+ " Intercept = " +Inter);
                //screen.println("Dist 1 = " +Dist1+ " Dist 2 = " +Dist2);
                //screen.println("xp1 = " +xp1+ "yp1 = " +yp1+ " xp2 = " +xp2+ " yp2 = " +yp2);
                //screen.println("x0 = " +x0+ "y0 = " +y0);
                //screen.println("delta in detector 1= "+d_det1);
                //screen.println("delta in detector 2= "+d_det2);
                //screen.println("delta in detector 3= "+d_det3);
                //screen.println(Q);
                //screen.println("Phi = " +Phi);
                //screen.println("x,y = " +x0 +", " +y0);
                //screen.println("BPhi = " +BPhi);
                //screen.println("Delta = " +d);
                //screen.println("Charge = " +q);
                //screen.println("LPhi = " +LPhi);
                //screen.println("x1= "+x1+ " y1= "+y1);
                //screen.println("x2= "+x2+ " y2= "+y2);
                //screen.println("x3= "+x3+ " y3= "+y3);
                //screen.println("Px= "+Prec [0][k]);
                //screen.println("Py= "+Prec [1][k]);
                //screen.println("Pz= "+Prec [2][k]);

                //static Histogram HiggsMass = new Histogram("Higgs Mass", 100, 0, 1000);

            }
            double InvMass = TrackMuonsBack.getInvMass(Prec);

            InvMass_array[i]=InvMass;

            screen.println("The invarient mass is " +InvMass);
            if (i < MuonData.size())
            {
                HiggsMass.fillh(InvMass);  // histogram the Higgs mass
            }

        }
        Import3.WriteToDisk(InvMass_array);
        HiggsMass.WriteToDisk(MuonData.size()); 
    }

    public static int getNumber()
    {
        return (Events*muons_per_event)+1;
    }

    public static void WriteToDisk(double []InvMass_array) throws IOException
    {
        String filename="..\\InvarientMass.csv"; //Creates a file with given name in directory one above the class directory
        FileWriter file1 = new FileWriter(filename); //this crates the file 
        PrintWriter outputFile = new PrintWriter (file1);

        outputFile.println("Event number , Inv Maass");
        //Double for loop to access all data in 3D array

        for (int i = 0; i<2022; i++)
        {
            //screen.println(j + " || " + trackOfMuon[i][0][j-1] + " || " + trackOfMuon[i][1][j-1]); // For debugging
            outputFile.println( i+1 + "," + InvMass_array [i]);
        }

        outputFile.close();
        screen.println("File written to disk in " +filename+ ".");
    }

}