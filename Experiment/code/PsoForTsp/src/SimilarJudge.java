
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


public class SimilarJudge {
	public static int Num=0;//Record the number of proteins after removal of impurities
	public static double[][] similar;
	public static String[] interactors;
	public static String[] go_iterm;
	
	
	 public List intersect(List ls, List ls2) {   
		 //the intersection operation
		          List list = new ArrayList(Arrays.asList(new Object[ls.size()]));   
		           Collections.copy(list, ls);   
		           list.retainAll(ls2);   //求二者交集，并将结果赋给list
		           return list;   
		       }   
		    
    public List union(List ls, List ls2) {   
    	//The union operation
		          List list = new ArrayList(Arrays.asList(new Object[ls.size()]));   
		          Collections.copy(list, ls);   
		          list.addAll(ls2);   
		           return list;   
		      } 
	
	/*
	 * 1.dataBuild is used to remove mismatched proteins (protein &amp;GO)
	 * liking:"Humaninteractors.txt", "HumanGO2.txt"
	 * 
	 */
	public  double[][] dataBuild(String interactorFile,String GoFile) throws IOException {
		
		
		BufferedReader br = new BufferedReader(new InputStreamReader(
		        new FileInputStream(interactorFile)));
	     String line=br.readLine();
	     //Interception of each protein, stored in an array
	     interactors =line.split(" "); 
		 
	     //Read GO data, where the data processed by the class of databuild
	     BufferedReader br2 = new BufferedReader(new InputStreamReader(
			        new FileInputStream(GoFile)));
	     String line2="";
	    
	     while ((line2 = br2.readLine())!= null) 
	     {//Remove the mismatched protein in the GO file, and then write to the file
			        	
	    	 go_iterm =line2.split(" ");
	    	 for (int i = 0; i < interactors.length; i++) {
				if (interactors[i].equals(go_iterm[0])) {
					 Num++;
					 fileWriter("Human3394Goentor.txt",interactors[i]+" ");//Production of non annotated proteins
					 fileWriter("Human3394_go.txt",line2+"\r\n");//Generate GO file
				}
			}
		}
	     br.close();
	     br2.close();
	     String goInputFile="Human3394_go.txt";
	     System.out.print(Num);//the number of protein
	     similar=new double[Num][Num];
	     BufferedReader br3 = new BufferedReader(new InputStreamReader(
			        new FileInputStream(goInputFile)));
	     BufferedReader br4 = new BufferedReader(new InputStreamReader(
			        new FileInputStream(goInputFile)));
	     //calculating the similarity matrix by comparing item of GO
	     for (int x = 0; x < Num; x++) 
	     {
	    	
	    	 String line3="";
	    	 String line4="";
	    	
	    	 if ((line3 = br3.readLine())!=null)
	    	 {
			for (int y = 0; y < Num; y++) 
		   {
				 
			if ((line4 = br4.readLine())!=null) 
			{
				
	 
			   if (y==x)
			   {
				   similar[x][y] = 1;
			     }
			   else
			   {
				   //Calculating the similarity
	        	   similar[y][x]=similar[x][y] = judge(line3, line4);
	        	  
			       
		       }
			}else {
			 br4 = new BufferedReader(new InputStreamReader(new FileInputStream(goInputFile)));
			 line4 = br4.readLine();
			}	
	   }
     }  	 
    }
     br3.close();
     br4.close();
		return similar;
}
	
	public double judge(String l1,String l2) {
		
			//Concrete steps to calculate similarity
		
			 String[] go_NUm1;
	  	     String[] go_NUm2;
  		     go_NUm1 =l1.split(" ");
			 go_NUm2 =l2.split(" ");
			 List list1=new ArrayList();
			 List list2=new ArrayList();
			 for (int i = 0; i < go_NUm1.length; i++) {
				
				 list1.add(go_NUm1[i]);
			}
	
			 for (int i = 0; i < go_NUm2.length; i++) {
				 
				 list2.add(go_NUm2[i]);
			}
			 List intersectList = intersect(list1, list2);   
			   
	         List unionList = union(list1, list2);   
	               
	           DecimalFormat df =new DecimalFormat("#.###");
	           double num1 = (double)intersectList.size();  //Get the number of intersection
	           double num2 = (double)unionList.size()-num1;  //Get the number of union
	           double similarNum=num1/num2;  //calculating the similarity
			   similarNum=Double.parseDouble(df.format(similarNum));
			   return similarNum;
		
		
	}
	
	//write to file
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
