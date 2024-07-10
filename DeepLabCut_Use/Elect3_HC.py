#cd C:\Users\Laura\Documents\DeepLabCut\DeepLabCut\conda-environments
# cd C:\Users\laura\DeepLabCut_Source\conda-environments
#conda activate dlc-windowsGPU
#conda activate DLC-GPU
# for desactivate
#conda deactivate
# for activate the ipython (perhaps don't needed with the IDLE
#ipython
import deeplabcut
import matplotlib
matplotlib.use('Agg') # for avoiding an error

# for creating a new project
config_path=deeplabcut.create_new_project('HomeCageElectrophy','Laura',List_video,working_directory='D:\Documents\Deeplabcut\Home_cage_Sleep_Stress',copy_videos=True)

# Nouveau projet sur Elect3
config_path= "D:\Documents\DeepLabCut_Root\HC_Elect3-Laura-2020-03-10\config.yaml"

config_path = 'D:\\DeepLabCut_Root\\HC_Elect3-Laura-2020-03-10\\config.yaml'
List_video=['D:\\DeepLabCut_Root\VideoTest\\200209_16_Hab1.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200209_21C_Hab1.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200213_16C_StressBL.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200218_17C_Stress_PS2.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200219_18C_Stress_PS1.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200223_19_Hab3.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200223_20_Hab3.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200224_19C_Stress_PS1.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200229_13_Hab2.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200229_14_Hab2.mp4',
            'D:\\DeepLabCut_Root\VideoTest\\200303_14C_Stress_PS2.mp4']

# for the white rats
List_video=['D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\Elect_BL_1session_121020_rat1.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\Elect_BL_1session_121020_rat2.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_BL1_22.10.20_rat2.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_BL2_23.10.20_rat3.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_BL3_06.11.20_rat2.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_JL1_09.11.20_rat3.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_JL1_26.10.20-BAD_rat2.mp4',
            'D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\LFP2_JL1_26.10.20-BAD_rat3.mp4'
            ]

# For adding new vids
deeplabcut.add_new_videos(config_path,List_video,copy_videos=True)

# for setting the config/settings (learning machine) you need several frame to set up manually, so that is extraction the frame
deeplabcut.extract_frames(config_path,'automatic','uniform',crop=False, userfeedback=True)

#enable easy labeling of all the extracted frames using an interactive GUI. The body parts to label (points of interest) should already have been named in the project’s configuration file (config.yaml, in Step 3).
deeplabcut.label_frames(config_path)

# checking the label
deeplabcut.check_labels(config_path)

# for training the network
deeplabcut.create_training_dataset(config_path,net_type='resnet_152',augmenter_type='imgaug', num_shuffles=1)

# used to train the network (very very long)
deeplabcut.train_network(config_path,shuffle=3)
deeplabcut.train_network(config_path)

# Evaluate the training network
deeplabcut.evaluate_network(config_path,plotting=True,Shuffles=[0,1,2,6,7,8])
deeplabcut.evaluate_network(config_path,plotting=True)


# Do a list of all vids for results
import os
path_vids = "D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\"

list_vids_full=[];
i = 1
for f in os.listdir(path_vids):
    name, ext = os.path.splitext(f)
    if ext == '.mp4':
        list_vids_full.append(path_vids+f)
        i = i+1
# Analyse videos and plotting data
deeplabcut.analyze_videos(config_path,['D:\Documents\Deeplabcut\Test_home_cage\Video2\HC1.mp4','D:\\DeepLabCut_Root\VideoTest\\200223_20_Hab3.mp4'],save_as_csv=True,videotype='.avi')
deeplabcut.analyze_videos(config_path,list_vids_full,save_as_csv=True,videotype='.mp4',shuffle=7)
deeplabcut.analyze_videos(config_path,list_vids_full,save_as_csv=True,videotype='.mp4')


# filter the predictions with a median filter (default; NEW as of 2.0.7+) or with a SARIMAX mode
deeplabcut.filterpredictions(config_path,['D:\\DeepLabCut_Root\VideoTest\\200213_16C_StressBL.mp4'])
deeplabcut.filterpredictions(config_path,list_vids_full)
deeplabcut.filterpredictions(config_path,list_vids_full,filtertype='median',shuffle=7)
deeplabcut.filterpredictions(config_path,list_vids_full,filtertype='arima',shuffle=7)

# This creates a new .h5 file with the ending _filtered that you can use in create_labeled_data and/or plot trajectories.

# The plotting components of this toolbox utilizes matplotlib therefore these plots can easily be customized by the end user. We also provide a function to plot the trajectory of the extracted poses across the analyzed video
deeplabcut.plot_trajectories(config_path,['D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\Elect_BL_1session_121020_rat1.mp4'])
deeplabcut.plot_trajectories(config_path,['D:\\DeepLabCut_Root\\AllVideoElect3_HC\\200209_16_Hab1.mp4'],shuffle=7,filtered = True)


# Plot the position on the video
deeplabcut.create_labeled_video(config_path,['D:\\DeepLabCut_Root\VideoTest\\200213_16C_StressBL.mp4'])

# We can add option :
# You can also optionally add a skeleton to connect points and/or add a history of points for visualization. To set the "trailing points" you need to pass trailpoints
# To draw a skeleton, you need to first define the pairs of connected nodes (in the config.yaml file) and set the skeleton color (in the config.yaml file). If you are using a project that was created before 2.0.7, you simply need to add these two items (skeleton and skeleton_color) to your config file (This addition is fully backwards compatible, so don't worry!).
# config skeleton
[['Nose', 'LeftEar'], ['Nose', 'RigthEar'], ['RigthEar', 'Hat'], ['LeftEar', 'Hat'], ['Hat','LeftShoulder'],
 ['Hat', 'RigthShoulder'], ['RigthShoulder', 'UpperBack'], ['LeftShoulder', 'UpperBack'], ['RigthShoulder', 'RigthHip'], ['LeftShoulder','LeftHip'],
 ['UpperBack', 'MiddleBack'], ['MiddleBack', 'LowerBack'], ['LeftHip', 'LowerBack'], ['RigthHip', 'LowerBack'], ['LowerBack','BaseTail'],
 ['BaseTail', 'MiddleTail'], ['MiddleTail', 'EndTail']]

- - Nose
  - Hat
  - LeftEar
  - RigthEar
- - LeftShoulder
  - RigthShoulder
  - UpperBack
- - MiddleBack
- - LeftHip
  - RigthHip
  - LowerBack
- - BaseTail
  - MiddleTail
  - EndTail

deeplabcut.create_labeled_video(config_path,['D:\\DeepLabCut_Root\VideoTest\\200223_20_Hab3.mp4'], filtered=True,trailpoints=10,draw_skeleton=True,save_frames=True)
deeplabcut.create_labeled_video(config_path,['D:\\DeepLabCut_Root\\AllVideoElect3_HC\\200209_16_Hab1.mp4'], filtered=True,trailpoints=5,draw_skeleton=True,shuffle=7)
deeplabcut.create_labeled_video(config_path,['D:\\Documents\\DeepLabCut_Root\\Video_rat_Blanc\\Elect_BL_1session_121020_rat1.mp4'], filtered=True)


# analyse skeleton ---- You can save the "skeleton" that was applied in create_labeled_videos for more computations. Namely, it extracts length and orientation of each "bone" of the skeleton as defined in the config.yaml file. You can use the function by:
deeplabcut.analyzeskeleton(config_path, list_vids_full, videotype='avi', shuffle=7, save_as_csv=True, destfolder=None)

# Refine the label
deeplabcut.refine_labels(config_path)

# merge the images after the refining
deeplabcut.merge_datasets(config_path)


# create a model comparison
# test with adding 3 suffles and changing the augmenter type
deeplabcut.create_training_model_comparison(config_path, num_shuffles = 1,
                                            net_types = ['resnet_50','resnet 101','resnet_152'],
                                            augmenter_types=['default','imgaug','deterministic'])

#########################################################################################
# GUI UPDATE
# they create a GUI for managing the project !!!
deeplabcut.launch_dlc()

# we can also continu to lauch the code as usual
# There are access to example with jupiter notebooks
jupyter notebook