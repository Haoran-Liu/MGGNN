cosine_value_threshold=0.05
patient_threshold=10
neighbors_threshold=4
max_percentage=0.05

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
control 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
sham 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
heme_0030 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
heme_0125 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
heme_0500 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
heme_1000 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage

Rscript /Users/haoran/Documents/MGGNN/data/sample_selection/Heme_spots_identification.r \
hemeHpx_1000 7 $cosine_value_threshold $patient_threshold $neighbors_threshold $max_percentage
