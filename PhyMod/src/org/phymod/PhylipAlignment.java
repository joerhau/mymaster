package org.phymod;

import org.phymod.io.Alignment;
import org.phymod.io.FileLoader;
import org.phymod.io.PartitionLoader;

public class PhylipAlignment extends Alignment {
	
	public PhylipAlignment(Species[] s) {
		spec = s;
	}
	
	public PhylipAlignment(String location, PartitionLoader part) {
		alignmentL = new FileLoader(location);
		alignmentL.open();
		this.partitionL = part;
		this.nrPartitions = part.nrPartitions;
		loadSpecies();
	}
	
	private void loadSpecies() {
		spec = new Species[alignmentL.l.size() - 1]; 
		for (int i = 0; i < alignmentL.l.size() - 1; i++) {
			String s = alignmentL.get(i + 1).replaceAll("\\s+", " "); 
			spec[i] = new Species(s.split(" ")[0].replaceAll("\\s+", ""), s.split(" ")[1].replaceAll("\\s+", ""), partitionL);
		}
	}

	
	public void toFile(String name) {
		FileLoader phylip = new FileLoader(name + ".phylip");
		FileLoader part = new FileLoader(name + ".part");
		phylip.write(this.phylipToString());
		part.write(this.partToString());
	}	
	
	public String phylipToString() {
		String s = " " + spec.length + " " + spec[0].length + "\n";
		for(int i = 0; i < spec.length; i++) {
			s += spec[i].name + " ";
			for(int j = 0; j < spec[i].partitions.length; j++) {
				s += spec[i].partitions[j].data;
			}
			s += "\n";
		}
		return s;
	}
	
	public String partToString() {
		String s = "";
		int start = 1;
		
		for(int j = 0; j < spec[0].partitions.length; j++) {
			int end = start + spec[0].partitions[j].data.length() - 1;
			s += spec[0].partitions[j].model + ", " + spec[0].partitions[j].name + " = " + start + "-" + end + "\n";
			start = end + 1;
		}
		return s;
	}
}
