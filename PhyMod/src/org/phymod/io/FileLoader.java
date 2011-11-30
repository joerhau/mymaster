package org.phymod.io;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.LinkedList;

public class FileLoader {
	public LinkedList<String> l;
	public File file;

	public FileLoader(String location) {
		this.file = new File(location);
	}
	
	public FileLoader(File f) {
		this.file = f;
	}

	public void write(String content) {
		try {
			if(file.exists()) {
				file.delete();
				file.createNewFile();
			} else file.createNewFile(); 
			
		    BufferedWriter out = new BufferedWriter(new FileWriter(file));
	    	out.write(content);
			out.close();
			
		} catch (Exception e) {
			System.out.println("Could nor open " + file.getAbsolutePath());
		}
	}

	public LinkedList<String> open() {
		try {
		    BufferedReader in = new BufferedReader(new FileReader(file));
		    l = new LinkedList<String>();
		    String tmp = "";
			while ((tmp = in.readLine()) != null) {
				if(tmp.trim().length() > 0) l.add(tmp);
			}
			in.close();
		} catch (Exception e) {
			System.out.println("Could not open " + file.getAbsolutePath());
		}
		
		return l;
	}
	
	public String line(int i) {
		return l.get(i);
	}
	
	public String toString() {
		String s = "";
		for(int i = 0; i < l.size(); i++) {
			s += l.get(i) + "\n";
		}
		return "File has " + l.size() + " nonempty lines\nContent is:\n" + s;
	}
}
