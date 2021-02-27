import java.util.Arrays;

public class Driver 
{
    public static void main( String[] args )
    {
        double y[] = {1000, 2000, 1000};
        double[][] v = { {1, 0,-1},{-1, 1,0},{1, -1,0},{-1, 0,-1}, {1, 0,-1},{-1, 1,0},{1, -1,0},{-1, 0,1}, {1,-1,0}  };
        double c[] = {5000.0, 50.0, 0.00005, 5.0, 5000.0, 50.0, 0.00005, 5.0, 0.05,};
    	
        Gillespie gillespie = new Gillespie(c, v, y, 3, 9, 0.0, 0.002);
        gillespie.run();
    }
}
