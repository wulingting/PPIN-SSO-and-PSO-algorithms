A total of three large folders are cleaned data,code,result.

/*****************************************************************************/
1.The file of cleand data:There are total 12 text files that all have been cleaned up,
corresponding to four species and each species includes the file of interaction, 
interactor and GO.
The data source is shown below:
---------------------------------------------------------------------------------
Abbreviation  |   Full name	                |           URL                  |
---------------------------------------------------------------------------------|
    DIP	      |Database of Interacting Proteins	|http://dip.doe-mbi.ucla.edu/dip/|
    GO	      |Gene Ontology Consortium	        |http://geneontology.org/        |
  UniProt     |Universal Protein	        |http://www.uniprot.org/         |
---------------------------------------------------------------------------------

/*****************************************************************************/

2.The file of code:This folder contains a detailed code of SSO and PSO algorithm.

2.1.The experimental environment for algorithm: operating system si 64 Windows8.1 
which is the professional edition; Inter processor (R) Pentium (R) CPU G3240 @ 3.10GHz;
running memory is 12GB; the data processing, the realization and running of the code are
main use of Java language in MyEclipse10.0.

2.2.The src file in the program contains databulid.java, FileWrite.java and so on.
Below I will give a detailed account of their general functions, of course, I also 
made the necessary comments in the code, so that readers can read easily.

(1)PSO.java:This is the main class of the project, the function includes the initialization 
   of the data, call the function of other classes, as well as a simple functional operation,
   and so on.Detailed steps are as follows:

   a.init():The first step is to read the data, and construct the similarity matrix and adjacency matrix;
   
   b.solve():The second step is to initial particle swarm,and then evolution;
  
   c.Optimization():The third step is to optimize the operation according to the similarity matrix and adjacency matrix;

   d.resultPrecess():The last step is to deal with the results, including filtering and merging the results, and calculate
     the value of the relevant issues, in order to evaluate the pros and cons of the algorithm.

(2)FileWriter.java:This class is a write operation to a file.

(3)getData.java:This class is to structure distance matrix.

(4)resultJudge.java:This class is to caculate the value of the relevant issues,inlucing Co,Se and so on.

(5)SimilarJudge.java:This class is to structure similarity matrix.

/*****************************************************************************/

2.The file of result:They are the results of PSO and SSO algorithm experiments, 
each algorithm corresponds to the results of four different species, as well as
the integration of the results of the four species.