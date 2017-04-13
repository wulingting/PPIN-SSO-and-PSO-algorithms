import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;


public class FilesWriter {

public void fileWriter(String fileName,String content){
		
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
