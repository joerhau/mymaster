package org.phymod;

public class Partition {
	public String model;
	public String data;
	public String name;
//	shall this partition be included in generated output
	public boolean masked;
	
	public Partition(String data, String name, String model) {
		this.model = model;
		this.name = name;
		this.data = data;
		this.masked = false;
	}
}
