package org.phymod;

import java.util.List;
import java.util.ArrayList;

public class Taxa {
	// name of this species
	public String name;
	// partitions list
//	private Partition[] partitions; 
	private List<Partition> partitions;
	// number of characters for this species
	public int length;
	// number of partitions available for this species
	public int nrPartitions;
	
	/**
	 * default constructor
	 * @param name name of the taxon
	 * @param partitions data assigned to this taxon
	 */
	public Taxa(String name, Partition[] partitions) {
		this.partitions = new ArrayList<Partition>(partitions.length);
		for(int i = 0; i < partitions.length; i++)
			this.partitions.add(partitions[i]);
		this.name = name;
		this.update();
	}
	
	/**
	 * update length, has to be called after anything within the partitions were updated
	 * @return
	 */
	private Taxa update() {
		int sum = 0;
		for (int i = 0; i < partitions.size(); i++)
			sum += partitions.get(i).data.length();
		this.length = sum;
		this.nrPartitions = partitions.size();
		
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
		this.partitions.remove(nr);
		return this.update();
	}
	
	public Partition getPartition(int nr) {
		return partitions.get(nr);
	}
	
	public String toString() {
		return name;
	}
}