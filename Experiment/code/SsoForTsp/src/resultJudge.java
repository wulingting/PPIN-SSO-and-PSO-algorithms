import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;


public class resultJudge {
	
	public void cohesion(double[][] similar,double[][] distence, String siteFile) throws IOException {
	    DecimalFormat df =new DecimalFormat("#.###");
		BufferedReader br1 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFile)));
		String line1="";
		int ModelNum=0;
		double CoSum=0.0;
		while ((line1=br1.readLine())!=null) {
			String[] site;
			double s,d;
			double Co=0.0;
			
			
			site=line1.split(",");
			int l=site.length;
			ModelNum++;
			for (int i = 0; i < l; i++) {
				
				if (i==l-1) {
					s=similar[Integer.parseInt(site[i])][Integer.parseInt(site[0])];
					d=distence[Integer.parseInt(site[i])][Integer.parseInt(site[0])];
				}else {
					s=similar[Integer.parseInt(site[i])][Integer.parseInt(site[i+1])];
					d=distence[Integer.parseInt(site[i])][Integer.parseInt(site[i+1])];
				}
			
				if (d!=0) {
					Co+=s+(1/d);    //Degree of polymerization formula, s is the similarity between proteins, D is the distance between proteins
				}else{
					Co+=s;
				}
			}
			CoSum +=(Co/l);//l is the number of proteins in the module
			
		}
		br1.close();
		System.out.print("The overall degree of polymerization is:"+df.format(CoSum/ModelNum)+"\r\n");//Total degree divided by module number
		
		
		
	}
	
	
	public void separation(double[][] similar,double[][] distence, String siteFile) throws IOException {
		/*
		 * In order to compare the two lines of the file, read the file two times, the Line1 and the line2 are always separated by one line;
		 * the Line3 records the first line for the last comparison
		 * 
		 */
		DecimalFormat df =new DecimalFormat("#.###");
		BufferedReader br1 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFile)));
		BufferedReader br2 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFile)));
		String line1="";
		String line2="";
		String line3=br2.readLine();
		int ModelNum1=0;
		int ModelNum2=0;
		double s,d;//Similarity and distance variables
		double SeSum=0.0;
		
		while ((line1=br1.readLine())!=null) {
			double Se=0.0;
			String[] site1;
			String[] site2;
			ModelNum1++;//Module number, module 1
			//System.out.print(line1+"\r\n");
			site1=line1.split(",");
			int l=site1.length;
			int m=0;
			if((line2=br2.readLine())!=null)
			{
				site2=line2.split(",");
				ModelNum2=ModelNum1+1;//module 2
				//System.out.print(line2+"\r\n");
			    
			for (int i = 0; i < l; i++) {
				for (int j = 0; j < site2.length; j++) {				
					s=similar[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];
					d=distence[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];	
					//System.out.print("s:"+s+"d:"+d);
					if (s!=0) {
						Se+= 1/s+d;  //The separation formula, s is similarity, D is the distance
					}else{
						Se+= d;
					}
					m++;
				}
	
				
			}	
			}
			else {
				site2=line3.split(",");
				for (int i = 0; i < l; i++) {
					for (int j = 0; j < site2.length; j++) {
						
						s=similar[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];
						d=distence[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];	
						//System.out.print("s:"+s+"d:"+d);
						if (s!=0) {
							Se+= 1/s+d;
						}else{
							Se+= d;
						}
						m++;
					}
		
				}	
				
			}
			if (ModelNum1==ModelNum2) {
				ModelNum2=1;
				//System.out.print("module"+ModelNum1+"and module"+ModelNum2+"\tSeparation: "+df.format(Se/m)+"\r\n");
				SeSum +=(Se/m);
			}else {
				//System.out.print("module"+ModelNum1+"and module"+ModelNum2+"\tSeparation: "+df.format(Se/m)+"\r\n");
				SeSum +=(Se/m);  //Total separation by module number
			}
			
			
		}
		
		br1.close();
		br2.close();
		System.out.print("The overall degree of separation between class:"+df.format(SeSum/ModelNum1)+"\r\n");
		
		
	}
	
	public String filterMerging(double[][] similar,double[][]distence,String siteFile) throws IOException {
		
		FilesWriter fl=new FilesWriter();
		String  siteFileString1 ="3394modelsitefile1.txt";
		String  siteFileString2 ="3394modelsitefile2.txt";
		
		DecimalFormat df =new DecimalFormat("#.###");
		BufferedReader br1 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFile)));
		BufferedReader br2 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFile)));
		String line1="";
		String line2="";
		String line3=br2.readLine();
		int ModelSize=0;//module size
		int modelNum=0;//the number of module
		double s;//Similarity and distance variables
		double SiSum=0.0;
		double Sijudge=0.0;
		double sijudgeSum=0.0;
		double disMark=0.005;//the degree of filtration 
		double SiMark=40.5;//the degree of combination
		int m=0;
		
		while ((line1=br1.readLine())!=null) {
		
			String[] site1;
			String[] site2;
			//System.out.print(line1+"\r\n");
			site1=line1.split(",");
			
			
			if((line2=br2.readLine())!=null)
			{
				site2=line2.split(",");
				ModelSize=site1.length<site2.length?site1.length:site2.length;
				//System.out.print(line2+"\r\n");
			    
			for (int i = 0; i < site1.length; i++) {
				for (int j = 0; j < site2.length; j++) {
					
					s=similar[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];
					//d=distence[Integer.parseInt(site1[i])][Integer.parseInt(site2[j])];	
					//System.out.print("s:"+s+"d:"+d);
					SiSum +=s;
					
				}	
			}
			Sijudge=SiSum/ModelSize;
			sijudgeSum+=Sijudge;
			m++;
			//System.out.print("\n Sijudge:"+Sijudge+"&"+m);
			
			
			}
			if (Sijudge>SiMark && line2!=null) {
				String mergeString=line1+line2+"\r\n";
				fl.fileWriter(siteFileString1, mergeString);//file record after merging
				line1=br1.readLine();
				line2=br2.readLine();
				
			} else {
				fl.fileWriter(siteFileString1, line1+"\r\n");
				//fl.fileWriter(siteFileString1, line2+"\r\n");

			}	
			modelNum++;
			
		}
		System.out.print("\nActual combining degree:"+df.format(sijudgeSum/modelNum));//实际总体合并度
		br1.close();
		br2.close();
		
		//Filtration process
	
		BufferedReader br3 = new BufferedReader(new InputStreamReader(
		        new FileInputStream(siteFileString1)));
		String linefiterString="";
		double d;
		double diSum=0.0;
		double disJudge=0.0;
		double disjudgeSum=0.0;
		int filtedNum=0;
		int modelNum2=0;
		int MSum=0;
		while ((linefiterString=br3.readLine())!=null) {
		   String[] site3;
		   site3=linefiterString.split(",");
		   modelNum2++;
		   for (int i = 0; i < site3.length; i++) {
				
			   if (i==site3.length-1) {
					
					d=distence[Integer.parseInt(site3[i])][Integer.parseInt(site3[0])];
				}else {
					
					d=distence[Integer.parseInt(site3[i])][Integer.parseInt(site3[i+1])];
				}
			   
			   diSum+=d;
				
			}	
		disJudge =site3.length/diSum ;//the number of node divided by length
		disjudgeSum+=disJudge;
		//System.out.print("\n dijudge"+disJudge);
		if (disJudge<disMark) {
			filtedNum+=site3.length;//If it is smaller than the density value, it is filtered and the number of filters is recorded
		}else {
			fl.fileWriter(siteFileString2, linefiterString+"\r\n");//If the filter value is greater than the value of the filtration, then write to the file
			MSum++;
		
		}
		
		
		
		}
		System.out.print("\nFinal module number:"+MSum);
		System.out.print("\nActual filtration density:"+df.format(disjudgeSum/modelNum2));
		System.out.println();
		return siteFileString2;
	}
	
}
