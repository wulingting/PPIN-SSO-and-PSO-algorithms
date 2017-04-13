import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import com.sun.org.apache.bcel.internal.generic.NEW;
import com.sun.org.apache.regexp.internal.recompile;


public class getData {
	
	
	public static int NUM=3394; //the number of proteins
	public static String[] interactors;
	public static String[] protein;
	public static int[][]  interaction=new int[NUM][NUM];//adjacency matrix
    public static double[][] dis=new double[NUM][NUM]; //distance matrix
	public static String[] Site;
	public static List intersectList;
	public static List unionList;
	public static double num1=0,num2=0,Distence=0;
	public static int lin_N=0;
	

	 public List intersect(List ls, List ls2) {   
		 //operation of intersection
		          List list = new ArrayList(Arrays.asList(new Object[ls.size()]));   
		           Collections.copy(list, ls);   
		           list.retainAll(ls2);   
		           return list;   
		       }   
		    
    public List union(List ls, List ls2) {   
    	//The union operation
		          List list = new ArrayList(Arrays.asList(new Object[ls.size()]));   
		          Collections.copy(list, ls);   
		          list.addAll(ls2);   
		           return list;   
		      } 
    
   //"Human3394Goentor.txt", "Humaninteractions.txt"
    public double[][] getDistence(String interactorFile,String interactionFile) throws IOException {
    	 BufferedReader br = new BufferedReader(new InputStreamReader(
			        new FileInputStream(interactorFile)));
		 String line="";
	     while ((line = br.readLine())!= null) 
	     {
	    	 //Interception of each protein, stored in an array       	
	    	 interactors =line.split(" "); //Protein species
		}
	     br.close();
	     BufferedReader str = new BufferedReader(new InputStreamReader(
			        new FileInputStream(interactionFile)));
		 while ((line = str.readLine())!= null) 
	     {
			 //Intercept interacting proteins, compare array of protein, and mark them in a two-dimensional array
			      protein =line.split(" ");
                  lin_N++;//Number of interactions
			      String Num="";
			       for (int i = 0; i < protein.length; i++) 
			    	 //Each interacting protein,liking:DIP-172E DIP-493N DIP-147N
			       {
					for (int j = 0; j < interactors.length; j++) 
						//Compared with one by one protein,DIP-98N DIP-96N DIP-966N DIP-95N DIP-953N DIP-952N
					{
							
			    	   if (protein[i].equals(interactors[j]))
			    	   {
						 
			    		   Num=Num+j+" ";//The order number of the interacting proteins including row and column
			    
			    	    }
					}
			    	   //System.out.println(protein[i]);
				  }
			      //Intercept coordinates and save in the array
			       
			       Site=Num.split(" ");
			      // adjacency matrix of the protein
			      if (Site.length<2) 
			      {
			    	  if(!Site[0].equals(""))
			    	  {
			    		//Diagonal assignment of adjacency matrix
			    		  interaction[Integer.parseInt(Site[0])][Integer.parseInt(Site[0])]=0;
			    		  //System.out.println(interaction[0][0]);
			          }
					
				  }else {
					 for (int i = 0; i < Site.length-1; i++) 
					 {
						 if(!Site[i].equals("")&& !Site[i+1].equals(""))
						 {
							 //Adjacency matrix element assignment
							 interaction[Integer.parseInt(Site[i])][Integer.parseInt(Site[i+1])]=1;
							 interaction[Integer.parseInt(Site[i+1])][Integer.parseInt(Site[i])]=1;
						 }

					 }
				  }
			      
		 }

		 str.close();
		 for (int i = 0; i < interaction.length; i++) 
		 {
			 List l1 = new ArrayList(); 
			 
			for (int n = 0; n < interaction.length; n++) 
			{
	
			 if (i==n) //Diagonal set to 1
			    {
					interaction[i][n]=1;
			    }
				if (interaction[i][n]==1) 
				{
					l1.add(n);//if the current row of proteins associated with each column of the protein, then set L1 1
				}
		    }
			
			for (int k = i+1; k < interaction.length; k++) 
			{
				 List l2 = new ArrayList(); 
				
			 for (int j = 0; j < interaction.length; j++) 
			 {

				if (k==j+1) 
				{
					interaction[k][j+1]=1;
					   
				}
				if (interaction[k][j]==1) 
				{
					l2.add(j);  //if the next line of protein is associated with a protein of  in each column, then  add to L2
				}
			 
				
			 }
			   intersectList =intersect(l1, l2);   	   
	            unionList =union(l1, l2);   
	   
	           DecimalFormat df =new DecimalFormat("#.###");
	           num1=(double)intersectList.size();  //Get the number of intersection
	           num2=(double)unionList.size()-num1;  //Get the number of Union
	           Distence=(num2-num1)/(num2+num1);  //Calculating distance
	           Distence= Double.parseDouble(df.format(Distence));
	           dis[i][k]=dis[k][i]=Distence;
			 
			}     
			 
		}
		return dis;
  
    	
	}
	
  //write to fille
  	public static void fileWriter(String fileName,String content){
  		
  	    File file = new File(fileName);
          BufferedWriter out=null;
  	    try {
  	    	out = new BufferedWriter(new FileWriter(file,true));
                 try {
   				out.append(content);
           	    out.flush();
                 } catch (IOException e) {
   				// TODO Auto-generated catch block
   				e.printStackTrace();
                 }
  	       
  	    } catch (IOException e1) {
  		// TODO Auto-generated catch block
  	    	e1.printStackTrace();
  	    }finally{
  	    	if(out!=null)
  				try {
  					out.close();
  				} catch (IOException e) {
  					// TODO Auto-generated catch block
  					e.printStackTrace();
  				}
  	    }
     	
  	}
	

}
