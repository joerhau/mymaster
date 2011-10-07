package org.phymod.io;

import org.phymod.Species;

public class Alignment {
	public Species[] spec;
	public int nrPartitions;
	
	public FileLoader alignmentL;
	public PartitionLoader partitionL;
	
	public String toString() {
		String s = "";
		for(int i = 0; i < spec.length; i++) {
			s += spec[i].name + "\n";
			for(int j = 0; j < spec[i].partitions.length; j++) {
				s += spec[i].partitions[j].data + "\n";
			}
		}
		return s;
	}
	
}
