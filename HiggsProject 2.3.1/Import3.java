import java.io.*;
import java.util.*;
public class Import3
{
    static PrintWriter screen = new PrintWriter( System.out, true);

    static String [] club = new String[20]; 
    static int Events = 2022; 
    static int entries_per_particle = 5;
    static int muons_per_event = 4;
    static double [] [] Finaldet1;
    static double [] [] Finaldet2;
    static double [] [] Finaldet3;

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

        Finaldet1 = new double [1+getNumber()] [3];
        Finaldet2 = new double [1+getNumber()] [3];
        Finaldet3 = new double [1+getNumber()] [3];

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
        Finaldet2 = TrackMuon.GetFinald2(); 
        Finaldet3 = TrackMuon.GetFinald3(); 

        for (int i = 0; i < MuonData.size(); i++) 
        {

            //coordinated of detector hits
            double x1= Finaldet1 [i][0];
            double y1= Finaldet1 [i][1];
            double z1= Finaldet1 [i][2];

            double x2= Finaldet2 [i][0];
            double y2= Finaldet2 [i][1];
            double z2= Finaldet2 [i][2];

            double x3= Finaldet3 [i][0];
            double y3= Finaldet3 [i][1];
            double z3= Finaldet3 [i][2];

            //double d = TrackMuonsBack.getd(x1, y1, z1, x2, y2, z2, x3, y3, z3);
            //double BPhi = TrackMuonsBack.getBPhi(x1, y1, z1);
            double B = 4; //B Field
            double R = 1.2;// Radius of the cavity

            //double LPhi = TrackMuonsBack.getLPhi(BPhi, d);

            //double [][]PT_array;
            //PT_array = new double [2022][4];
            //PT_array[i][0] = TrackMuonsBack.getP_T(d, B, R);
            //PT_array[i][1] = TrackMuonsBack.getP_T(d, B, R);
            //PT_array[i][2] = TrackMuonsBack.getP_T(d, B, R);
            //PT_array[i][3] = TrackMuonsBack.getP_T(d, B, R);

            double q= MuonStore[i][4];//carge of the muon
            double LPhi = TrackMuonsBack.getLPhi(x1, y1);
            double Pperp1 = TrackMuonsBack.getPperp(x1, y1);
            //double Pperp2 = TrackMuonsBack.getPperp(x2, y2);
            //double Pperp3 = TrackMuonsBack.getPperp(x3, y3);
            double d_det1 = TrackMuonsBack.getdelta(Pperp1, q, R, B);
            //double d_det2 = TrackMuonsBack.getdelta(Pperp2, q, R, B);
            //double d_det3 = TrackMuonsBack.getdelta(Pperp3, q, R, B);
            double BPhi = TrackMuonsBack.getPhi(LPhi, d_det1);

            /*if (d_det1 > 0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }
            if (d_det2 > 0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }
            if (d_det3 > 0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }
            if (d_det1 < -0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }
            if (d_det2 < -0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }
            if (d_det3 < -0.006)
            {
                screen.println("Low momentum muon, abandon it.");
                continue;
            }*/

            int muonsPerEvent = 4;

            double [][]PT1_array;
            PT1_array = new double [ MuonData.size()][muonsPerEvent];
            //double [][]PT2_array;
            //PT2_array = new double [ MuonData.size()][muonsPerEvent];
            //double [][]PT3_array;
            //PT3_array = new double [ MuonData.size()][muonsPerEvent];

            double [][]Pz_array;
            Pz_array = new double [ MuonData.size()][muonsPerEvent];

            double [][]E_array;
            E_array = new double [ MuonData.size()][muonsPerEvent];

            for (int j = 0; j < muonsPerEvent; j++) 
            {
                PT1_array[i][j] = TrackMuonsBack.getPperp(x1, y1);
                //PT2_array[i][j] = TrackMuonsBack.getPperp(x2, y2);
                //PT3_array[i][j] = TrackMuonsBack.getPperp(x3, y3);

                Pz_array[i][j] = MuonStore[i][j];

                E_array[i][j] = MuonStore[i][j];

                //screen.println("Pz= " +Pz_array[i][j]+ ", PT 1= " +PT1_array[i][j]);//+ ", PT 2= " +PT2_array[i][j]+ ", PT 3= " +PT3_array[i][j]);
            }

            /*double [][]Q_array;
            Q_array = new double [2022][4];
            Q_array[i][0] = TrackMuonsBack.getQ(d);
            Q_array[i][1] = TrackMuonsBack.getQ(d);
            Q_array[i][2] = TrackMuonsBack.getQ(d);
            Q_array[i][3] = TrackMuonsBack.getQ(d);*/

            //double parentMass = TrackMuonsBack.getOriginalMass(PT_array[i][0],PT_array[i][1],PT_array[i][2],PT_array[i][3],
            //Pz_array[i][0],Pz_array[i][1],Pz_array[i][2],Pz_array[i][3],
            //E_array[i][0],E_array[i][1],E_array[i][2],E_array[i][3]);

            //double parentAbsMom = TrackMuonsBack.getAbsMomentum(PT_array[i][0],PT_array[i][1],PT_array[i][2],PT_array[i][3],
            //Pz_array[i][0],Pz_array[i][1],Pz_array[i][2],Pz_array[i][3]);

            //double parentCharge = TrackMuonsBack.getOriginalQ(Q_array[i][0],Q_array[i][1],Q_array[i][2],Q_array[i][3]);

            //screen.println("The parent particle's mass is " +parentMass+ ", the parent particle's absolute momentum is " 
            //+parentAbsMom+ ".");
            //Q = TrackMuonsBack.getQ(d_det1);

            //screen.println("delta in detector 1= "+d_det1);
            //screen.println("delta in detector 2= "+d_det2);
            //screen.println("delta in detector 3= "+d_det3);
            //screen.println(Q);
            //screen.println(LPhi);
            //screen.println(BPhi);
            screen.println("x1= "+x1+ "y1= "+y1+ "z1= "+z1);
            screen.println("x2= "+x2+ "y2= "+y2+ "z2= "+z2);
            screen.println("x3= "+x3+ "y3= "+y3+ "z3= "+z3);
        }
    }

    public static int getNumber()
    {
        return Events*muons_per_event+1;
    }

}