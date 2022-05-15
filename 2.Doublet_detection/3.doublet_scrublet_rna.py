import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import sys
import pandas as pd
import pynndescent
import datetime
import os

starting_time = datetime.datetime.now()

print("Starting job at " + starting_time.strftime("%Y-%m-%d %H:%M:%S"))

data_path=os.path.join(os.getcwd(),"Documents/multiome_tonsil_Lucia/2.doublet_detection/tmp/each_sample")
print(data_path)

#data_path = "/Users/odelacalle/Documents/multiome_tonsil_Lucia/2.doublet_detection/tmp/"
#os.chdir(data_path)

dirlist = [ item for item in os.listdir(data_path) if os.path.isdir(os.path.join(data_path, item)) ]
print ("directory list:", dirlist)


for filename in os.listdir(data_path):
    print(filename)
    f = os.path.join(data_path, filename)
    print(f)
    # checking if it is a file
    if os.path.isdir(f):
      
        counts_matrix = scipy.io.mmread(os.path.join(f, "rna_matrix.mtx")).T
        
        folder_figures = os.path.join(f,"Figures")
        pred_doublet=os.path.join(f,"scrublet_text")
        data_frames= (os.path.join(f,"data_frames"))
        if not os.path.exists(folder_figures):
            os.mkdir(folder_figures)
        else:
            pass
        if not os.path.exists(pred_doublet):
            os.mkdir(pred_doublet)
        else:
            pass

        if not os.path.exists(data_frames):
            os.mkdir(data_frames)
        else:
            pass

        path_to_barcodes=os.path.join(f,"barcodes.tsv")
        barcodes=pd.read_csv(path_to_barcodes,header=None)
        scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.05)

        print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))


        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 75, n_prin_comps = 30)
        doublet_scores = np.round(doublet_scores, decimals = 3)

        scrub.plot_histogram()[0].savefig(os.path.join(folder_figures,"scrublet_histogram.png"))

        print('Running UMAP...')

        scrub.set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))
        scrub.plot_embedding('UMAP', order_points=True)[0].savefig(os.path.join(folder_figures,"UMAP_doublets.png"))


        # Write out predicted doublet scores 
        np.savetxt(fname = os.path.join(pred_doublet,"scrublet_doublet_scores.txt"), X = doublet_scores)
        np.savetxt(fname = os.path.join(pred_doublet,"scrublet_predicted_doublets.txt"), X = predicted_doublets)

        print("Creating data frame...")

        print(doublet_scores)
        print(predicted_doublets)

        #arr_barcodes = pd.DataFrame.to_numpy(barcodes)
        
        df_doublet_scores = pd.DataFrame(doublet_scores,columns=["doublet_scores"])
        df_predicted_doublets=pd.DataFrame(predicted_doublets,columns=["predicted_doublets"])

        
        print(df_doublet_scores)
        print(df_predicted_doublets)
        print(barcodes)
        print(type(df_doublet_scores))
        print(type(df_predicted_doublets))


       

        #scrublet_doubl_dict = {"barcodes": barcodes, "scrublet_doublet_scores": df_doublet_scores, "scrublet_predicted_doublet": df_predicted_doublets}
        #print(scrublet_doubl_dict)
        frame=[barcodes,df_doublet_scores,df_predicted_doublets]
        scrublet_doubl_dict1 =pd.concat(frame,axis=1)
        scrublet_doubl_dict1_new= scrublet_doubl_dict1.rename(columns={0:'barcodes'})
    
        print("****** printing dataframe ****")
        print(type(scrublet_doubl_dict1))
        print(scrublet_doubl_dict1)
        print("****** printing NEW dataframe ****")
        print(scrublet_doubl_dict1_new)



        df_predicted_doublets.to_csv(os.path.join(data_frames,"scrublet_prediction.csv"),index = False)
        df_doublet_scores.to_csv(os.path.join(data_frames,"scrublet_doublet_score.csv"), index = False)
        scrublet_doubl_dict1_new.to_csv(os.path.join(data_frames,"scrublet_doublet_prediction.csv"), index = False)

ending_time = datetime.datetime.now()



print("Ending job at " + ending_time.strftime("%Y-%m-%d %H:%M:%S"))
print('Done.')

# Run
#print("Running scrublet...")


#expected_doublet_rate = 0.056
# for key in all_data.keys():
#     print(key)
#     # Scrublet
#     scrub = scr.Scrublet(all_data[key][0], expected_doublet_rate = 0.6)
#     if "atac" in key:
#         doublet_scores, predicted_doublets = scrub.scrub_doublets(log_transform = True, min_counts = 2, min_cells = 3, min_gene_variability_pctl = 70, n_prin_comps = 50)
#     elif "rna" in key:
#         doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts = 2, min_cells = 3, min_gene_variability_pctl = 75, n_prin_comps = 30)
#     doublet_scores = np.round(doublet_scores, decimals = 3)
    
    
#     # Create DataFrame
#     scrublet_doubl_dict = {"barcodes": all_data[key][1][0].values, "scrublet_doublet_scores": doublet_scores, "scrublet_predicted_doublet": predicted_doublets}
#     scrublet_doubl_df = pd.DataFrame(scrublet_doubl_dict)
#     all_data[key].append(scrublet_doubl_df)
    
    
#     # Plot
#     scrub.set_embedding("UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist = 0.3))
#     hists = scrub.plot_histogram()
#     umaps = scrub.plot_embedding("UMAP", order_points = True)
    

#     # Save
#     scrublet_doubl_df.to_csv("/Users/odelacalle/Documents/multiome_tonsil_Lucia/results/scrublet/scrublet_doublet_prediction_{}.csv".format(key), index = False)
#     hists[0].savefig("/Users/odelacalle/Documents/multiome_tonsil_Lucia/2.doublet_detection/tmp/histograms/scrublet_doublet_prediction_histograms_{}.png".format(key), dpi = 100)
#     umaps[0].savefig("/Users/odelacalle/Documents/multiome_tonsil_Lucia/2.doublet_detection/tmp/umaps/scrublet_doublet_prediction_umaps_{}".format(key), dpi = 100)
