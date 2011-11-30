package org.phymod;

public class Taxa {
	public String name;
	public Partition[] partitions;
	public int length;
	
	/**
	 * default constructor
	 * @param name name of the taxon
	 * @param partitions data assigned to this taxon
	 */
	public Taxa(String name, Partition[] partitions) {
		this.partitions = partitions;
		this.name = name;
		this.update();
	}
	
	/**
	 * update length, has to be called after anything within the partitions were updated
	 * @return
	 */
	private Taxa update() {
		int sum = 0;
		for (int i = 0; i < partitions.length; i++)
			sum += partitions[i].data.length();
		this.length = sum;
		return this;
	}
	
	/** 
	 * adds a partition to an already existing taxon
	 * @param part
	 * @return
	 */
	public Taxa addPartition(Partition part) {
		Partition[] ps = new Partition[this.partitions.length + 1];
		
		for(int i=0; i < this.partitions.length; i++)
			ps[i] = this.partitions[i];
		ps[this.partitions.length] = part;
		this.partitions = ps;
		
		return this.update();
	}
	
	public String toString() {
		return name;
	}
}