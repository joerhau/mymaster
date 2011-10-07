package org.phymod;

public class Partition {
	public String model;
	public String data;
	public String name;

	
	public Partition(String data, String name, int start, int end, String model) {
		this.model = model;
		this.name = name + "." + start + "-" + end;
		this.data = data;
		
	}
	
	public Partition(String data, String name, String model) {
		this.model = model;
		this.name = name;
		this.data = data;
	}
}
