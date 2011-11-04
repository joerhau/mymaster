package org.phymod;

import org.phymod.io.PartitionLoader;

public class Species {
	public String name;
	public Partition[] partitions;
	public int length;
	
	public Species(String name, int length) {
		this.name = name;
		this.partitions = new Partition[length];
		
	}
	
	public Species(String string, String string2, PartitionLoader part) {
		this.name = string;
		initPartitions(part, string2);
	}

	public void initPartitions(PartitionLoader p, String data) {
		this.partitions = new Partition[p.nrPartitions];
		this.length = data.length();
		
		for(int i = 0; i < p.nrPartitions; i++) {
			partitions[i] = new Partition(data.substring(p.start[i] - 1, p.end[i]), p.names[i], p.start[i], p.end[i], p.models[i]);
		}
	}
	
	public void replace(Species s, int source, int dest) {
		partitions[dest] = new Partition(s.partitions[source].data, s.partitions[source].name, s.partitions[source].model);
		length += s.partitions[source].data.length();
	}
	
	public void addPartition(Partition part) {
		Partition[] ps = new Partition[this.partitions.length + 1];
		for(int i=0; i < this.partitions.length; i++)
			ps[i] = this.partitions[i];
		ps[this.partitions.length] = part;
		this.partitions = ps;
		this.length++;
	}
	
	public String toString() {
		return name;
	}
}
