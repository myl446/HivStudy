package com.cn;
/*
 *This class is mainly used to calculate the distance matrix of species, input: DNA sequence file.
 *@author:Ma Yuanlin
 *@2019-7-20
 */
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ComputerdistanceMatrix {

	public static void main(String[] args) throws IOException {
		ComputerdistanceMatrix csf=new ComputerdistanceMatrix();
		//String FILE_IN = "F:\\mayuanlinstudy\\HIVdatastudy\\datacg1625";
		String FILE_IN = "F:\\mayuanlinstudy\\HIVdatastudy\\PureSequences825";
		// String FILE_IN = "F:\\mayuanlinstudy\\HIVdatastudy\\5596data";

		for(int k=8;k<9;k++){
			 File file = new File("F:\\mayuanlinstudy\\HIVdatastudy\\PureSequences825"+"\\ManhK="+k+".txt");
				if (!file.exists()) {
					file.createNewFile();
				} // if file doesnt exists, then create it
		    FileWriter fileWritter = new FileWriter(file.getAbsoluteFile());
			BufferedWriter bufferWritter = new BufferedWriter(fileWritter);//
			
		 List<String> list = new ArrayList<String>();
			File f = new File(FILE_IN);
			list = getFileList(f); 
			
			
			/*for (String l : list) {
				bufferWritter.write(l.substring((l.lastIndexOf("\\"))+1,l.lastIndexOf("."))+"\t\t");
			}*/
			
			// bufferWritter.write("\t"+list.size()+"\n");
			long startTime=System.currentTimeMillis();
			Map feature=new LinkedHashMap<String, double[]>();
			
		    for (String l : list) {/
		    	   String tempString = null;
		    	   double[] aa=null;
		    	  // System.out.println(l);
					//System.out.println(list.size());
					File rfile=new File(l);
					//bufferWritter.write(l.substring((l.lastIndexOf("\\"))+1,l.lastIndexOf("."))+"\t\t");
					String speName=l.substring((l.lastIndexOf("\\"))+1,l.lastIndexOf("."));
					//System.out.println(l);
					BufferedReader br = new BufferedReader(new FileReader(rfile));
					int fileLong=0;
					String temps="";
					while ((tempString = br.readLine())!= null) {
						if (tempString.startsWith(">")) {
							tempString=br.readLine();
							
					} 
				
						temps=temps+tempString;
					fileLong=fileLong+tempString.length();
					}
					
					br.close();
					aa=csf.computWeightKermasend(temps, k);
					feature.put(speName, aa);
	
				}
		    
		    Set<String> keySet=feature.keySet();
		    for(String name: keySet){
		    	bufferWritter.write(name+"\t");
		    	double[] bb=(double[])feature.get(name);
		    	Set<String> keySet1=feature.keySet();
		    	 for(String name1: keySet){
		    		 double[] cc=(double[])feature.get(name1);
		    		 double distance=csf.computManHdistance(bb, cc);
		    		 bufferWritter.write(distance+"\t");
		    		 }
		    	 bufferWritter.write("\n");
		    	
		    }					
			 
		    long endTime=System.currentTimeMillis();
		    System.out.println("k="+k+"   矩阵运行时间"+(endTime-startTime)/1000+"秒");
		    bufferWritter.flush();
		    bufferWritter.close();
	}	
}
	
	
	
	double computManHdistance(double[] a,double[] b){
		double sum=0.0;
		for(int i=0;i<a.length;i++){
			sum=sum+Math.abs(a[i]-b[i]);
			}
		//System.out.println("sum="+sum);
		return sum;
		}
	
	public double[]  computWeightKermasend(String sequence, int kermas){
		
		
		int[] strnum=this.seqTonumvec(sequence);
		double[] feater=new double[(int)Math.pow(4, kermas)];
		List<Map> list=new ArrayList<>();
	    for(int i=0;i<strnum.length-kermas+1;i++){
	    	int m=0;
	    	for(int j=0;j<kermas;j++){
	    		m=4*m+strnum[i+j];
	          }
	    	Map map=new HashMap();
	    	map.put(m,i+1);
	    	list.add(map);
	    }
	     Map temp=mapCombine(list);
	             Iterator<Map.Entry<Integer, List>> it = temp.entrySet().iterator();
	             while (it.hasNext()) {
	                  Map.Entry<Integer, List> entry = it.next();
	                  List templist=entry.getValue();
	                  int[] temparray=new int[templist.size()];
	                  double sum=0.0;
	                  //***************************************************************************************************
	                
	                  for(int k=0;k<templist.size();k++){
	                	  temparray[k]=((Integer)templist.get(k)).intValue();
	                	  sum=sum+temparray[k];
	                	   }
	                  feater[entry.getKey()]=1.0*sum/((strnum.length-kermas+1)*strnum.length);//
	                  
	             }
	             temp.clear();
	             temp=null;
	             return feater;
		
	}
	public static int[] seqTonumvec(String str){
		int[] numseq=new int[str.length()];
		char ch=' ';
		for(int i=0;i<str.length();i++){
			ch=str.charAt(i);
			switch (ch) {
			case 'A':
				numseq[i]=0;
				break;
			case 'C':
				numseq[i]=1;
				break;
			case 'G':
				numseq[i]=2;
				break;
			case 'T':
				numseq[i]=3;
				break;

			default:
				numseq[i]=(int)(Math.random()*4);
				break;
			}
		}
		return numseq;
		
	}
	
	public int numofNoAGCT(String str){
		int[] numseq=new int[str.length()];
		int n=0;
		char ch=' ';
		for(int i=0;i<str.length();i++){
			ch=str.charAt(i);
			switch (ch) {
			case 'A':
				numseq[i]=0;
				break;
			case 'C':
				numseq[i]=1;
				break;
			case 'G':
				numseq[i]=2;
				break;
			case 'T':
				numseq[i]=3;
				break;

			default:
				numseq[i]=(int)(Math.random()*4);
				n++;
				break;
			}
		}
		return n;
		
	}
	
	
	public static Map mapCombine(List<Map> list){
		Map<Object,List> map=new HashMap<>();
		
		for(Map m:list){
			
			Iterator<Object> it=m.keySet().iterator();
			while(it.hasNext()){
				Object key=it.next();
				if(!map.containsKey(key)){
					List newList=new ArrayList<>();
					newList.add(m.get(key));
					map.put(key, newList);
					newList=null;
					
				}else{
					map.get(key).add(m.get(key));
				}
				
			}
			m.clear();
			m=null;
		}
		return map;
		
		
	}

	public static List<String> getFileList(File file) { 
	       List<String> result = new ArrayList<String>();  
  if (!file.isDirectory()) {  
      System.out.println(file.getAbsolutePath());  
      result.add(file.getAbsolutePath());  
    } else {  
      File[] directoryList = file.listFiles(new FileFilter() {  
         public boolean accept(File file) {  
             if (file.isFile() && file.getName().indexOf("ffn")> -1) {  
                 return true;  
             } else {  
                 return false;  
             }  
          }  
     }
      );  
  for (int i = 0; i < directoryList.length; i++) {  
         result.add(directoryList[i].getPath());  
         }  
           }  
     return result;  
 } 
	public static String getFileNameNoEx(String filename){
		if(filename!=null&&(filename.length()>0)){
			int dot=filename.lastIndexOf('.');
			int cha=filename.lastIndexOf('\\');
			if((dot>-1)&&(dot<(filename.length()))){
				return filename.substring(cha+1, dot);
			}
			
		}
		return filename;
	}

	
	  public static String format4(double value) {
		 return new Formatter().format("%.6f", value).toString();
		}
	  
	  public static List getFilefoldList(String filepath){
		  File file = new File(filepath);
		  File[] fileList = file.listFiles();
		  List wjjList = new ArrayList<File>();
		  for (int i = 0; i < fileList.length; i++) {
		     if (fileList[i].isDirectory()) {
		          wjjList .add(fileList[i]);
		     }
		  }
		  return wjjList;
		  
	  }
}
