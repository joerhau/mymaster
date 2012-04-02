package org.phymod;

import java.util.List;
import java.util.ArrayList;

public class Taxa {
	// name of this species
	public String name;
	// partitions list
//	private Partition[] partitions; 
	public List<Partition> partitions;
//	number of characters for this species
	public int length;
//	number of partitions available for this species
	public int nrPartitions;
//	shall this species be included in generated output
	public boolean masked;
	
	/**
	 * default constructor
	 * @param name name of the taxon
	 * @param partitions data assigned to this taxon
	 */
	public Taxa(String name, Partition[] partitions) {
		this.partitions = new ArrayList<Partition>(partitions.length);
		for(int i = 0; i < partitions.length; i++)
			this.partitions.add(partitions[i]);
		this.masked = false;
		this.name = name;
		this.update();
	}
	
	/**
	 * update length, has to be called after anything within the partitions were updated
	 * @return
	 */
	public Taxa update() {
		nrPartitions = 0;
		length = 0;
		for (int i = 0; i < partitions.size(); i++) {
			nrPartitions += !partitions.get(i).masked ? 1 : 0;
			length += !partitions.get(i).masked ? partitions.get(i).data.length() : 0;
		}
		return this;
	}
	
	/** 
	 * adds a partition to an already existing taxon
	 * @param part
	 * @return
	 */
	public Taxa addPartition(Partition part) {
		this.partitions.add(part);
		return this.update();
	}
	
	public Taxa removePartition(int nr) {
		partitions.get(nr).masked = true;
		return update();
	}
	
	public Partition getPartition(int nr) {
		return partitions.get(nr);
	}
	
	public String toString() {
		return name;
	}
}