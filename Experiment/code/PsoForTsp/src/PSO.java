/**
*-------------------------------------------------------------------------------------------------------------
* PsoForTsp algorithm by Wu Lingting        
* Programmed by Zheng Xianghan, Chen Riqing, Ye Shaozhen at College of Mathematics and Computer Science, Fuzhou University, Fuzhou, China
* Last revised: April 2017     
*--------------------------------------------------------------------------------------------------------------
* Paper -- Citation Details: 
* 1）Xianghan Zheng, Lingting Wu1, Shaozhen Ye1, Riqing Chen 
* Simplified Swarm Optimization based Function Modules Detection in Protein-to-Protein Interaction Networks
**/	
import java.io.BufferedReader;  
import java.io.FileInputStream;  
import java.io.IOException;  
import java.io.InputStreamReader;  
import java.text.DecimalFormat;
import java.util.ArrayList;  
import java.util.List;
import java.util.Random;  
import org.omg.CORBA.PUBLIC_MEMBER;
import com.sun.org.apache.bcel.internal.generic.NEW;
	  
	public class PSO {  
		getData getdata= new getData();
		SimilarJudge similarJudge=new SimilarJudge();
		DecimalFormat df =new DecimalFormat("#.###");
	    
	    private double w;  //weight
	    private int MAX_GEN;// The number of iterations  
	    private int scale;// Population size or population scale
	    private int proteinNum; // Number of proteins,that is coding length 
	    
	    private int bestNum;  
	    private int t;// current t generation    
	  
	    private double[][] distance; // Distance matrix  
	    private double[][] similar;//Similar matrix
	      
	    private double[][] oPopulation;// particle swarm   
	    private ArrayList<ArrayList<SO>> listV;// The initial exchange sequence of every division particles  
	  
	    private double[][] Pd;// The generation of the best solution for a particle  
	    private double[] vPd;// Evaluation value of solution  
	  
	    private double[] Pgd;// The best solution of the whole particle swarm, each particle can remember the best solution  
	    private double vPgd;// The evaluation value of the best solution 
	    private int bestT;// The bestT generation of the best solution    
	  
	    private double[] fitness;// The fitness of the population, indicating the fitness of each individual in the population
	    private double c1=2,c2;
	     
	    private double[] record;
	    private double fitSum=0.0;//The length of the path
	    private double gama=0.65;//The cutoff value of parameter
	    private double cutoffDis;//The cutoff value of distance
	    private double cutoffSim;//The cutoff value of similarity
	    
	    private int pro_M;//the number of hypothetical module
	    private double Ds=0.07;//The density of the module
	    private int lin_N;//The number of interaction
	  
	    private List list1;
	    private Random random;  
	  
	    public PSO() {  
	  
	    }  
	  
	    /**
	     *  
	     * @param n 
	     *            The number of protein 
	     * @param g 
	     *            The number of operation 
	     * @param s
	     *            Population size
	     * @param w 
	     *            weight 
	     **/  
	    public PSO(int n, int g, int s) {  
	        this.proteinNum = n;  
	        this.MAX_GEN = g;  
	        this.scale = s; 
	        
	        
	    }  
	  
	    @SuppressWarnings("resource")  
	    /** 
	     * The initialization of PSO algorithm 
	     * @throws IOException 
	     */  
	    private void init() throws IOException {  
	        // Read data 
	        distance = new double[proteinNum][proteinNum];  
	        similar = new double[proteinNum][proteinNum];  
	        //distance=getdata.getDistence("interactors.txt", "interactions.txt");

	        similar=similarJudge.dataBuild("Humaninteractors.txt", "HumanGO2.txt");//Original tor file + original GO file
	        //similar=similarJudge.dataBuild("Scereinteractors.txt","ScereGo2.txt" );
	        distance=getdata.getDistence("Human3394Goentor.txt", "Humaninteractions.txt");
	       // distance=getdata.getDistence("Scere3394Goentor.txt", "Scereinteractions.txt");//tor file after removing impurity+ original tion file
	        System.out.println("Start2....");
	        
	        oPopulation = new double[scale][proteinNum];  
	        fitness = new double[scale];
	        Pd = new double[scale][proteinNum];  // The best solution for a particle in each generation
	        vPd = new double[scale];  //Evaluation value of solution
	        Pgd = new double[proteinNum];  //The best solution of the whole particle swarm, each particle can remember the best solution 
	        vPgd = Integer.MAX_VALUE; // The evaluation value of the best solution of the whole particle swarm
	        bestT = 0;  
	        t = 0;  
	  
	        random = new Random(System.currentTimeMillis());  
	     
	    }  
	  
	    // Initialize the population to generate the unique number
	    void initGroup() {  
	        int i, j, k;  
	        for (k = 0; k < scale; k++)// Population size 
	        {  
	            oPopulation[k][0] = random.nextInt(65535) % proteinNum;  
	            for (i = 1; i < proteinNum;)// Spatial dimension, here is the number of proteins
	            {  
	                oPopulation[k][i] = random.nextInt(65535) % proteinNum;  
	                for (j = 0; j < i; j++) {  
	                    if (oPopulation[k][i] == oPopulation[k][j]) {  
	                        break;  
	                    }  
	                }  
	                if (j == i) {  
	                    i++;  
	                }  
	            }  
	        }             
	    }  
	  
	    void initListV() {  
	        int ra;  
	        int raA;  
	        int raB;  
	  
	        listV = new ArrayList<ArrayList<SO>>();  
	  
	        for (int i = 0; i < scale; i++) {  
	            ArrayList<SO> list = new ArrayList<SO>();  
	            ra = random.nextInt(65535) % proteinNum;  //Generate the number between 0 and proteinNum
	            for (int j = 0; j < ra; j++) {  
	                raA = random.nextInt(65535) % proteinNum;  
	                raB = random.nextInt(65535) % proteinNum;  
	                while (raA == raB) {  
	                    raB = random.nextInt(65535) % proteinNum;  
	                }  
	  
	                // raA and raB are differernt  
	                SO s = new SO(raA, raB);  
	                list.add(s);  
	            }  
	  
	            listV.add(list);  
	        }  
	    }  
	  
	    public double evaluate(double[] chr) 
	    {  
	        // 0123  
	        double len = 0;  
	        
	        /*Here is the TSP algorithm based on the traveling salesman,
	        that is, coding, the starting city, the city 1, the city of the city of 2 ...n 
	        (where the city is equivalent to the protein)*/ 	         
	        for (int i = 1; i < proteinNum; i++) {  
	            len += distance[(int)chr[i - 1]][(int)chr[i]];  
	        }  
	        // City n, starting city  
	        len += distance[(int)chr[proteinNum - 1]][(int)chr[0]];  
	        return len;  
	    }  
	  
	    // coding after a basic sequence on arr exchange coding  
	    public void add(double[] arr, ArrayList<SO> list) {  
	        double temp = -1;  
	        SO s;  
	        for (int i = 0; i < list.size(); i++) {  
	            s = list.get(i);  
	            temp = arr[s.getX()];  
	            arr[s.getX()] = arr[s.getY()];  
	            arr[s.getY()] = temp;  
	        }  
	    }  
	  
	    // the basic exchange sequence of two codes, such as A-B=SS
	    public ArrayList<SO> minus(double[] a, double[] b) {  
	        double[] temp = b.clone();  
	        int index;  
	        // Commutator 
	        SO s;  
	        // Exchange sequence
	        ArrayList<SO> list = new ArrayList<SO>();  
	        for (int i = 0; i < proteinNum; i++) {  
	            if (a[i] != temp[i]) {  
	                // Find the same value as a[i] index in the array of temp  
	                index = findNum(temp, a[i]);  
	                // exchange the index of i and index value in the arrary of temp 
	                changeIndex(temp, i, index);  
	                // Remember the commutator 
	                s = new SO(i, index);  
	                // save the commutator
	                list.add(s);  
	            }  
	        }  
	        return list;  
	    }  
	  
	    //Find num in the array of arr and return the corresponding index value 
	    public int findNum(double[] arr, double num) {  
	        int index = -1;  
	        for (int i = 0; i < proteinNum; i++) {  
	            if (arr[i] == num) {  
	                index = i;  
	                break;  
	            }  
	        }  
	        return index;  
	    }  
	  
	    // exchange the index of index1 and index2 value in the array of arr  
	    public void changeIndex(double[] arr, int index1, int index2) {  
	        double temp = arr[index1];  
	        arr[index1] = arr[index2];  
	        arr[index2] = temp;  
	    }  
	  
	    // Two-dimensional array copy 
	    public void copyarray(double[][] from, double[][] to) {  
	        for (int i = 0; i < scale; i++) {  
	            for (int j = 0; j < proteinNum; j++) {  
	                to[i][j] = from[i][j];  
	            }  
	        }  
	    }  
	  
	    // One-dimensional array copy
	    public void copyarrayNum(double[] from, double[] to) {  
	        for (int i = 0; i < proteinNum; i++) {  
	            to[i] = from[i];  
	        }  
	    }  
	      
	    public void evolution() {  
	        int i, j, k;  
	        int len = 0;  
	        float ra = 0f;  
	  
	        ArrayList<SO> Vi;  
	          
	        // One iteration    
	        for (t = 0; t < MAX_GEN; t++) {  
	            // For each particle
	            for (i = 0; i < scale; i++) {  
	                 if(i==bestNum) continue;  
	                ArrayList<SO> Vii = new ArrayList<SO>();  	                
	                // Updating speed  	                
	                Vi = listV.get(i);  
	                
	                w=0.9-0.4*t/MAX_GEN;
	                // the size of Vi * w (rounding) exchange sequence 
	                len = (int) (Vi.size() * w);  
	                
	                for (j = 0; j < len; j++) {  
	                    Vii.add(Vi.get(j));  
	                }  
	  
	                // Pid-Xid  
	                ArrayList<SO> a = minus(Pd[i], oPopulation[i]);  
	                ra = random.nextFloat();  
	  
	                // ra(Pid-Xid)+  
	                len = (int)(a.size()*ra);  
	                
	                for (j = 0; j < len; j++) {  
	                    Vii.add(a.get(j));  
	                }  
	  
	                // Pid-Xid  
	                ArrayList<SO> b = minus(Pgd, oPopulation[i]);  
	                //ra = random.nextFloat();  
	               
	                c2=t/MAX_GEN;
	                // ra(Pid-Xid)+  
	                len = (int)(b.size()*c2);   
	                for (j = 0; j < len; j++) {  
	                    SO tt= b.get(j);  
	                    Vii.add(tt);  
	                }  
	  
	                // save the updating Vii  
	                listV.add(i, Vii);  
	  
	                // Updating location  
	                // Xid¡¯=Xid+Vid  
	                add(oPopulation[i], Vii);  
	            }  
	  
	            //  Calculating the fitness of the new particle swarm and selecting the best solution
	            for (k = 0; k < scale; k++) 
	            {  
	                fitness[k] = evaluate(oPopulation[k]);  
	                if (vPd[k] > fitness[k]) 
	                {  
	                    vPd[k] = fitness[k];  
	                    copyarrayNum(oPopulation[k], Pd[k]);  
	                    bestNum=k;  
	                }  
	                if (vPgd > vPd[k]) 
	                {  
	                    System.out.println("Optimum length"+Double.parseDouble(df.format(vPgd))+" generation£º"+bestT);  
	                    bestT = t;  
	                    vPgd = vPd[k];  
	                  
	                    copyarrayNum(Pd[k], Pgd);  
	                }  
	            }         
	        }  
	    }  
	  
	    public void solve() {  
	        int i;  
	        int k;  
	  
	        initGroup();  
	        initListV(); //Each particle is a vector of D dimensions, and each dimension is a two-dimensional list 
	  
	        // Each particle remembers its best solution 
	        copyarray(oPopulation, Pd);  //OPopulation is a randomly generated particle swarm
	        // Calculating the fitness of the new particle swarm and selecting the best solution,Fitness[max]
	        for (k = 0; k < scale; k++) {  
	            fitness[k] = evaluate(oPopulation[k]); //referring to the problem of TSP, calculating the distance len from the initial city to traverse all cities and back to the original location 
	            vPd[k] = fitness[k];  
	            if (vPgd > vPd[k]) {  
	                vPgd = vPd[k];  
	                copyarrayNum(Pd[k], Pgd);  
	                bestNum=k;  //Find the best solution corresponding to the population number bestNum
	            }  
	        }  
	  
	        // output  
	        System.out.println("Initial particle swarm...");  
	        for (k = 0; k < scale; k++) {  
	            for (i = 0; i < proteinNum; i++) {  
	                System.out.print((int)oPopulation[k][i] + ",");  
	            }  
	            System.out.println();  
	            System.out.println("----" + Double.parseDouble(df.format(fitness[k])));  
	        }  
	  
	        // evolution 
	        evolution();  
	  
	        // output
	        System.out.println("final particle swarm...");  
	        for (k = 0; k < scale; k++) {  
	            for (i = 0; i < proteinNum; i++) {  
	                System.out.print((int)oPopulation[k][i] + ",");  
	            }  
	            System.out.println();  
	            System.out.println("----KKKK" + Double.parseDouble(df.format(fitness[k])));  
	            fitSum += fitness[k];//total path length
	          
	        }  
	          
	        System.out.println("The optimal path generation:");  
	        System.out.println(bestT);  
	        System.out.println("The length of optimal path");  
	        System.out.println(Double.parseDouble(df.format(vPgd)));  
	        System.out.println("optimal path:");  
	        for (i = 0; i < proteinNum; i++) {  
	            System.out.print((int)Pgd[i] + ",");  
	        }  
	        }  
	    
	    public  void Optimization()
	    {
			
	    	double judgeDis;
	    	double judgeSim;
	    	double judge;
	    	double e;
	    	
	        list1=new ArrayList();
	    	cutoffDis=(fitSum/scale)/proteinNum;//The cutoff value = average path length * parameter
	    	
	    	for (int i = 0; i < proteinNum; i++)
	    	{
				
	    		if (i==proteinNum-1)
	    		{
	    			judgeDis=distance[(int)Pgd[0]][(int)Pgd[proteinNum-1]];//The distance between nodes
	    			judgeSim=similar[(int)Pgd[0]][(int)Pgd[proteinNum-1]];//Get the similarity between nodes
	    			judge=judgeDis*judgeSim;
	    		}
	    		else {
				    judgeDis=distance[(int)Pgd[i]][(int)Pgd[i+1]];//The distance between nodes
				    judgeSim=similar[(int)Pgd[i]][(int)Pgd[i+1]];//Get the similarity between nodes
				    judge=judgeDis*judgeSim;
				    }
	    		   System.out.print("*"+judge);

	    			if (judge>=0.08) 
	    		{
	    			//judgeSum=0.0;
					list1.add(i);//The path of getting truncated coordinates
				}
			}
	    	/*System.out.println("Truncated coordinates£º ");
	    	for (int i = 0; i < list1.size(); i++) 
	    	{
	    		System.out.println(" "+Integer.parseInt(String.valueOf(list1.get(i))));
			}*/
	    	
	    	
	    	
			if (list1.isEmpty())
			{
				for (int i = 0; i < proteinNum; i++) 
				{
					System.out.print((int)Pgd[i] + ",");
				}
				
			}else {
				int temp2=0;
				for (int i = 0; i < list1.size(); i++) 
				{
					
					int temp=Integer.parseInt(String.valueOf(list1.get(i)));
					
				    if (i==0) 
				    {
				    	System.out.print("\n"+"module"+(i+1)+": ");
				    	for (int j = 0; j <= temp ; j++) 
				    	{
				    		System.out.print((int)Pgd[j] + ",");
				    		
						}
				    	temp2=temp+1;
					}else 
                         {
						
							System.out.print("\n"+"module"+(i+1)+": ");
							for (int j = temp2; j <= temp ; j++) 
							{
								System.out.print((int)Pgd[j] + ",");
					        }
							temp2=temp+1;
					     }
				    if (i==list1.size()-1) 
							    {
									System.out.print("\n"+"module"+(list1.size()+1)+": ");
									
									for (int j = temp2; j <proteinNum; j++) 
									{
							    		System.out.print((int)Pgd[j] + ",");
									}
									
									
								}
				       
					}
				
				}
			
			
			
		}
	
	    
	    
	    FilesWriter fr=new FilesWriter();
	    //the manipulation function of module coordinate
	    public String  ModelSite() {
	    	String siteFileNameString="modelSiteFile3394.txt";
	    	
	    	int filteredNum=0;
	    	int filterSize=2;//Granularity of filtration
	    
	    	if (list1.isEmpty())
			{
	    		String lineString1 ="";
				for (int i = 0; i < proteinNum; i++) 
				{
				
					lineString1=(int)Pgd[i]+",";
					
				}
				
				fr.fileWriter(siteFileNameString,lineString1+"\r\n");
			}else {
				

				int temp2=0;
				for (int i = 0; i < list1.size(); i++) 
				{
					
					int temp=Integer.parseInt(String.valueOf(list1.get(i)));
					
				    if (i==0) 
				    {
				    	String lineString2 ="";
				    	int mark1=0;//Number of proteins in the module
				    	for (int j = 0; j <= temp ; j++) 
				    	{
				    	
				    		lineString2+=(int)Pgd[j]+",";//Get coordinate of protein
				    		mark1++;
						}
				  if (mark1<filterSize) {//If the value is lower than the threshold filtering, filtering
					filteredNum+=mark1;
					lineString2 =null;
					temp2=temp+1;//If the value is lower than the threshold filtering, filtering
				}else {
				  
				    	fr.fileWriter(siteFileNameString,lineString2+"\r\n");//write to the file by line
				    	temp2=temp+1;
				    }
					}else 
                         {
						    String lineString3 ="";
						    int mark2=0;//Number of proteins in the module
							for (int j = temp2; j <= temp ; j++) 
							{
						
								lineString3+=(int)Pgd[j]+",";
								mark2++;
					        }
					
							if (mark2<filterSize) {
								filteredNum+=mark2;
							    lineString3 =null;
							    temp2=temp+1;
								
							}else {
							  
							    	fr.fileWriter(siteFileNameString,lineString3+"\r\n");//write to the file by line
							    	temp2=temp+1;
                         }
					     }
				    if (i==list1.size()-1) 
							    {
							
				    	            String lineString4 ="";
				    	            int mark3=0;//Number of proteins in the module
									for (int j = temp2; j <proteinNum; j++) 
									{
							    	
							    		lineString4+=(int)Pgd[j]+",";
										mark3++;
									}
						
									if (mark3<filterSize) {
										filteredNum+=mark3;
										lineString4 =null;
										temp2=temp+1;
									}else {
									  
									    	fr.fileWriter(siteFileNameString,lineString4+"\r\n");//write to the file by line
									  
							    }
									
								}
			
					   }

				}

	    	System.out.print("\n The number of filtered proteins£º"+ filteredNum);
	    	
	    	
	    	
	    	
			return siteFileNameString;
		}
	    
	    
	    public void resultPrecess() throws Exception {
	    	
	    	String file=ModelSite();
	    	resultJudge rj=new resultJudge();
	    	String siteProcessedFile=rj.filterMerging(similar, distance, file);
	    	rj.cohesion(similar, distance, siteProcessedFile); //computing the value of the degree of polymerization in the module(Co) 
	    	rj.separation(similar, distance, siteProcessedFile); //computing the value of the deviation degree between modules 
			
		}
	    
	    
	    
	    public class SO {  
	    	  
	        private int x;  
	        private int y;  
	          
	        public SO(int x,int y)  
	        {  
	            this.x=x;  
	            this.y=y;  
	        }  
	          
	        public int getX() {  
	            return x;  
	        }  
	        public void setX(int x) {  
	            this.x = x;  
	        }  
	        public int getY() {  
	            return y;  
	        }  
	        public void setY(int y) {  
	            this.y = y;  
	        }  
	          
	        public void print()  
	        {  
	            System.out.println("x:"+this.x+" y:"+this.y);  
	        }  
	    }  
	  
	    /** 
	     * @param args 
	     * @throws Exception 
	     */  
	    public static void main(String[] args) throws Exception {  
	        System.out.println("Start....");  
	        long  startTime=System.currentTimeMillis();	      
	         PSO pso = new PSO(1447, 500, 100);//number of protein, number of iterations, population size,Mouse
	        //PSO pso = new PSO(3394, 500, 100);//number of protein, number of iterations, population size,Human 
	        //PSO pso = new PSO(269, 500, 100);//number of protein, number of iterations, population size,Fruitfly
	        // PSO pso = new PSO(2325, 500, 100);//number of protein, number of iterations, population size,Scere13
	        pso.init();  
	        System.out.println("start1 finish");
	        pso.solve();
	        System.out.println("start2 finish");
	        pso.Optimization();
	        System.out.println("start3 finish");
	        pso.resultPrecess();	     
	        System.out.println("start4 finish");
	        long  endTime=System.currentTimeMillis();
	        System.out.println("\nTime= "+(endTime-startTime)/1000);  //running time of the algorithm
	    }  
	}  

