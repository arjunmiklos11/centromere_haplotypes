for chr in chr{1..22};
do
	echo "module reset; module load miniconda; conda activate /home/sl2749/.conda/envs/tensorqtl; python input_generator.py $chr" >> run.input_generator.sh
done	
