package org.phymod.io;

public class PartitionLoader extends FileLoader{
	public int nrPartitions;
	public String[] names;
	public String[] models;
	public int[] start;
	public int[] end;
	
	public PartitionLoader(String location) {
		super(location);
		nrPartitions = l.size();
		parse();
	}
	
	private void parse() {
		names = new String[nrPartitions];
		models = new String[nrPartitions];
		start = new int[nrPartitions];
		end = new int[nrPartitions];
		
		for (int i = 0; i < nrPartitions; i++) {
			String str = l.get(i).replaceAll("\\s+", "");
			models[i] = str.split(",")[0];
			names[i] = str.split(",")[1].split("=")[0];
			start[i] = Integer.parseInt(str.split(",")[1].split("=")[1].split("-")[0]);
			end[i] = Integer.parseInt(str.split(",")[1].split("=")[1].split("-")[1]);
		}
	}
	
	public String toString() {
		String s = "";
		for(int i = 0; i < nrPartitions; i++) {
			s += models[i] + ", " + names[i] + " = " + start[i] + "-" + end[i] + "\n";
		}
		return "Folling partitions have ben parsed: \n" + s;
	}
}
