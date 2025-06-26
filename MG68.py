import multiprocessing, os, sys, time, imagej, shutil, pyautogui, keyboard
import scipy.stats as stats, pandas as pd, numpy as np, cv2 as cv
import cellprofiler_core.utilities.java
import cellprofiler_core.preferences
import cellprofiler_core.pipeline
from pywinauto.application import Application
from multiprocessing import Process
from scipy.stats import ttest_ind
from scyjava import jimport
#KaranMalhotra_HebertLab_6/15/24

np.seterr(all = "ignore")
#ignore all irrelevant warnings...

flaggedimages, goodimages = [], []
screenWidth, screenHeight = pyautogui.size()
#globalvariables

def inputfunction(anotheriter):
    # this function gets input data from the user, and stores/passes it back to "main"
    # it is relatively simple, and easy to follow. Just be careful of changing the order of the return variables, as they are passed to downstream functions based on their order.

    print("\nBefore you run this program, make sure your folders are clear...errors will be raised if other information is in the folders you input...")

    if not anotheriter:

        # check whether this is "anotheriter" (i.e. a call to this function beyond the first original call). "anotheriter" would be true if you edited images and are running them through the pipeline again.

        global inputfiles
        global outputfolder
        global correcthresh
        global integralint
        global filename
        global inputfolder
        global yeskey
        global booldefaultanalysis
        global densitybool
        global microglialdataframename
        global Pipepath
        # store variables
        # I am lazy, so I defined many of these "input" variables as global, so you can reference them in other functions. In practice, this is usually highly discouraged. For our purposes, since this isn't a super long program (and hence, variables aren't likely to be confused), it shouldn't matter, but if you want to improve the code, you can use OOP to implement classes/attributes.

        print("\n Choose analysis...")
        print("\n-> Sholl/CellP Analysis ('s')...")

        densityinput = input("-> Density Analysis ('d')...\n")

        if densityinput in ['s', 'yes', 'sholl', 'S', 'Sholl', 'c', 'CellP']:

            densitybool = False

        else:

            densitybool = True

        # check what type of analysis the user wants (either sholl/cellp or density)

        while True:

            try:

                inputfolder = input("\nEnter the path of your input folder...\n")
                inputfolder = inputfolder[1:len(inputfolder) - 1]
                inputfiles = os.scandir(inputfolder)
                # be careful with "os.scandir" if you implement it elsewhere. It gets file paths, but it is an iterator (you can only iterate over its items once, before it is exhausted).
                # make a list from os.scandir, and then make copies, or subset that list. Do not iterate over the same os.scandir object.
                break

            except:

                print("\nYour path is invalid: enter with quotes, and ensure you have the full path!...\n")

        # get path of input folder (try and except is a way to catch an error in the path name. If there is an error in copy/paste (and the path is invalid), the user can attempt again.

        keycheck = input("\nDo your files have a keyword? ('y' or 'n')...\n")
        if keycheck in ['y', 'yes', 'yah', 'yeah', 'yup']:

            keyword = input("\ntype the keyword, exactly as it appears...\n")
            inputfiles = [entry.path for entry in inputfiles if keyword in entry.path]
            yeskey = [True, keyword]

        else:

            inputfiles = [entry.path for entry in inputfiles]
            yeskey = [False, False]

        # check for keywords, and if so, filter the input files for them + store note of the key for future reference (in another run)

        if len(inputfiles) == 0:

            print("\nerror, no files found, quitting!\n")
            quit()

        # if no files, then quit/error

        outputfolder = input("\nEnter the path for output folder to save in...\n")
        outputfolder = outputfolder[1:len(outputfolder) - 1]
        # get the output folder

        filename = input("\nEnter name to save the files as...\n")
        # store the output file names. This name will be appended to the file names

        Pipepath = input("\nEnter the path of your pipeline...\n")
        # store the path of the pipeline. Unfortunately, there is no way to check validity (and catch errors), until the pipeline actually is running.

        if not densitybool:

            # if we aren't doing density analysis, then take paramters for the cellp/sholl analysis.
            # if you decide to add parameters you want to modulate on each run, in addition to adding them to the cellprofiler function, add them down here.
            # if they are numbers (i.e. a threshold value, make sure you convert them to float/integer. "input" as a function returns a string.

            correcthresh = (float(input("\nEnter the correctionthreshold value...\n")))
            integralint = int(input("\nEnter the integralint value...\n"))
            skipfiles = input("\nDo you want to skip file editing?... ('s' or 'e')\n")
            # take in correctionthreshold, integralvalue for filtering, and book for skipping editing.

            if skipfiles in ['s', 'skip', 'yes', 'y']:

                print("\nskipping editing...\n")
                yeskey.append(True)

            else:

                print("\nediting on, pls watch for output...\n")
                yeskey.append(False)

            # decide whether to skip editing. In practice, unless there is a surely known artifact, or you want to define a region for analysis, I don't think the editing feature should be used.

            dataanalyzersettingmarker = input("Do you want to use default analysis settings, or set your own ('d' or 'o')\n")
            if dataanalyzersettingmarker in ['d', 'default', 'yes', 'y']:

                booldefaultanalysis = True

            else:

                booldefaultanalysis = False

            # I added functionality for not using default analysis settings (which is implemented in the dataorganizer processing, and essentially cross-checks cellprofiler and imageJ, filterinng out cells that aren't overlapped enough between the two programs, based on nuclear location/size). By default this setting is turned off (commented out), because its too unreliable (cellp and imagej have differing identifcation methods, so its hard to get true overlap). If you want to make this setting relevant, you will need to change it in "dataorganizer."

            microglialdataframename = filename + "DataSheet" + "_ii" + str(integralint) + "_ct_" + str(correcthresh) + ".xlsx"
            # set filename here, so it can be referenced by other functions

        else:

            # if we are doing density analysis

            correcthresh = 0
            integralint = 0
            # need to set these to something, since returning them (could also do it at the beginning of the function)
            booldefaultanalysis = False
            yeskey.append(False)
            microglialdataframename = filename + "DensityDataSheet.xlsx"
            # set filename here, so it can be referenced by other functions

    else:

        # if we are in a run besides the original instance. (if "anotheriter" is true)

        inputfiles = [entry.path for entry in os.scandir(inputfolder)]
        filename = input("\nEnter name to save the files as\n")
        keyword = input("\ntype the keyword you used, exactly as it appears...\n")
        inputfiles = [entry for entry in inputfiles if keyword in entry]
        yeskey = [True, keyword]
        # access files if they want to retain settings, only need to re-ask for file names (with the keyword they used for saving editing images), so we aren't overwriting the new files with new settings.

        if len(inputfiles) == 0:

            print("\nno files found, quitting!\n")
            quit()

        # again quit if files are invalid

        skipfiles = input("\nDo you want to skip file editing?... ('s' or 'e')\n")
        if skipfiles in ['s', 'skip', 'yes', 'y']:

            print("\nskipping editing...\n")
            yeskey.append(True)

        else:

            print("\nediting on, pls watch for output...\n")
            yeskey.append(False)

        # same as above

        microglialdataframename = filename + "DataSheet" + "_ii" + str(integralint) + "_ct_" + str(correcthresh) + ".xlsx"
        # set name

    return [inputfiles, inputfolder, Pipepath, outputfolder, correcthresh, integralint, filename, yeskey, densitybool]
    #return saved variables

def cellprofiler(inputfilelist, Pipepath, outputfolder, correcthresh, integralint, filename, densitybool):

    # run cellprofiler

    currenttime = time.time()
    # record current time

    cellprofiler_core.preferences.set_headless()
    cellprofiler_core.utilities.java.start_java()
    # start up cell profiler + JVM

    Pipeline = cellprofiler_core.pipeline.Pipeline()
    Pipeline.load(Pipepath[1:len(Pipepath) - 1])
    # load the pipeline

    print("\nrunning analysis...\n")

    if densitybool:

        # if doing density analysis

        Pipeline.modules()[8].setting(10).set_value("Elsewhere...|" + outputfolder)
        # setting output (the numbers here are dependent on your cellprofiler pipeline). I've included info in the guide on how they work + how to check what settings correspond to what number.

        Pipeline.read_file_list(inputfilelist)
        output_measurements = Pipeline.run()
        # read the file list + run the pipeline

        MicrogliaCenterXoutput = output_measurements.get_all_measurements("PossibleNuclei", "Location_Center_X")
        MicrogliaCenterYoutput = output_measurements.get_all_measurements("PossibleNuclei", "Location_Center_Y")
        # get locations and store (returns an iterable of all measurements).

        cellprofiler_core.utilities.java.stop_java()
        # kill the jvm

        MicrogliaX = []
        MicrogliaY = []
        namelist = []
        namecounter = 0
        # storage

        for centerx, centery in zip(MicrogliaCenterXoutput, MicrogliaCenterYoutput):

            i = 1

            while i <= len(centerx):

                namelist.append(inputfilelist[namecounter])
                i += 1

            MicrogliaX.extend(centerx)
            MicrogliaY.extend(centery)
            namecounter += 1

        # append numbers + names to lists

        MicroglialDataFrame = {"FileName": namelist, "MicrogliaXPosition": MicrogliaX, "MicrogliaYPosition": MicrogliaY}
        # create a dictionary of filenames + measurements

        microglialdataframename = outputfolder + "\\" + filename + "DensityDataSheet.xlsx"
        pd.DataFrame(MicroglialDataFrame).to_excel(microglialdataframename, index_label="Index")
        # make a pandas dataframe + write it to excel, with index values for later filtering uses.

        print("cycle complete   |   " + str(round(((time.time() - currenttime) / 60), 3)) + " mins")
        print("\nstarting tracing...\n")

    else:

        # if doing sholl/cellp

        Pipeline.modules()[7].setting(16).set_value(integralint)
        # integratedintensity
        Pipeline.modules()[8].setting(15).set_value(correcthresh)
        # correctthreshold
        Pipeline.modules()[20].setting(10).set_value("Elsewhere...|" + outputfolder)
        # set image output names

        Pipeline.read_file_list(inputfilelist)
        output_measurements = Pipeline.run()
        # run pipeline

        MicrogliaCD68 = output_measurements.get_all_measurements("TrueMicroglia", "Intensity_MeanIntensity_OrigRed")
        MicrogliaCD68Integrated = output_measurements.get_all_measurements("TrueMicroglia", "Intensity_IntegratedIntensity_OrigRed")
        MicrogliaArea = output_measurements.get_all_measurements("TrueMicroglia", "AreaShape_Area")
        MicrogliaMajorAxis = output_measurements.get_all_measurements("TrueMicroglia", "AreaShape_MajorAxisLength")
        MicrogliaCenterX = output_measurements.get_all_measurements("TrueMicroglia", "Location_Center_X")
        MicrogliaCenterY = output_measurements.get_all_measurements("TrueMicroglia", "Location_Center_Y")
        NucleiArea = output_measurements.get_all_measurements("TrueNuclei", "AreaShape_Area")
        NucleiCD68 = output_measurements.get_all_measurements("TrueNuclei", "Intensity_MeanIntensity_OrigRed")
        NucleiCD68Integrated = output_measurements.get_all_measurements("TrueNuclei", "Intensity_IntegratedIntensity_OrigRed")
        NucleiCenterX = output_measurements.get_all_measurements("TrueNuclei", "Location_Center_X")
        NucleiCenterY = output_measurements.get_all_measurements("TrueNuclei", "Location_Center_Y")
        # get relevant measurements

        cellprofiler_core.utilities.java.stop_java()
        # kill the jvm

        MicrogliaAreaList = []
        MicrogliaCD68List = []
        NucleiAreaList = []
        NucleiCD68List = []
        MicrogliaMajorAxisList = []
        MicrogliaCD68IntegratedList = []
        NucleiCD68IntegratedList = []
        MicrogliaCenterXList = []
        MicrogliaCenterYList = []
        NucleiCenterXList = []
        NucleiCenterYList = []
        namelist = []
        # create lists

        namecounter = 0
        for marea, mcd68, narea, ncd68, mmajoraxis, mcd68int, ncd68int, mcenterx, mcentery, ncenterx, ncentery in zip(
                MicrogliaArea, MicrogliaCD68, NucleiArea, NucleiCD68, MicrogliaMajorAxis, MicrogliaCD68Integrated,
                NucleiCD68Integrated, MicrogliaCenterX, MicrogliaCenterY, NucleiCenterX, NucleiCenterY):

            i = 1

            while i <= len(marea):

                namelist.append(inputfilelist[namecounter])
                i += 1

            MicrogliaAreaList.extend(marea)
            MicrogliaCD68List.extend(mcd68)
            NucleiAreaList.extend(narea)
            NucleiCD68List.extend(ncd68)
            MicrogliaMajorAxisList.extend(mmajoraxis)
            MicrogliaCD68IntegratedList.extend(mcd68int)
            NucleiCD68IntegratedList.extend(ncd68int)
            MicrogliaCenterXList.extend(mcenterx)
            MicrogliaCenterYList.extend(mcentery)
            NucleiCenterXList.extend(ncenterx)
            NucleiCenterYList.extend(ncentery)
            namecounter += 1
        # name counter is unncessary--just for my own checks

        nc_over_mc_area = [x / y for x, y in zip(NucleiAreaList, MicrogliaAreaList)]
        maxis_area = [x / y for x, y in zip(MicrogliaMajorAxisList, MicrogliaAreaList)]
        # for division measurements

        MicroglialDataFrame = {"FileName": namelist, "MicrogliaCD68MeanIntensity": MicrogliaCD68List,
                               "MicrogliaCD68IntegratedIntensity": MicrogliaCD68IntegratedList,
                               "MicrogliaArea": MicrogliaAreaList,
                               "MicrogliaMajorAxisLength": MicrogliaMajorAxisList,
                               "MajorAxis/MicrogliaArea": maxis_area,
                               "NucleiArea": NucleiAreaList,
                               "NucleiCD68MeanIntensity": NucleiCD68List,
                               "NucleiCD68IntegratedIntensity": NucleiCD68IntegratedList,
                               "NucleiArea/MicrogliaArea": nc_over_mc_area,
                               "NucleiXPosition": NucleiCenterXList,
                               "NucleiYPosition": NucleiCenterYList,
                               "MicrogliaXPosition": MicrogliaCenterXList,
                               "MicrogliaYPosition": MicrogliaCenterYList}
        # dictionary for raw values

        MicroglialAverages = {"FileName": namelist,
                              "MicrogliaCD68MeanIntensity": [np.average(x) for x in MicrogliaCD68List],
                              "MicrogliaCD68IntegratedIntensity": [np.average(x) for x in MicrogliaCD68IntegratedList],
                              "MicrogliaArea": [np.average(x) for x in MicrogliaAreaList],
                              "MicrogliaMajorAxisLength": [np.average(x) for x in MicrogliaMajorAxisList],
                              "MajorAxis/MicrogliaArea": [np.average(x) for x in maxis_area],
                              "NucleiArea": [np.average(x) for x in NucleiAreaList],
                              "NucleiCD68MeanIntensity": [np.average(x) for x in NucleiCD68List],
                              "NucleiCD68IntegratedIntensity": [np.average(x) for x in NucleiCD68IntegratedList],
                              "NucleiArea/MicrogliaArea": [np.average(x) for x in nc_over_mc_area]}
        # dictionary for averages across all images

        MicroglialDataFrame = pd.DataFrame(MicroglialDataFrame)
        MicroglialAverages = pd.DataFrame(MicroglialAverages)
        # create dataframes

        microglialdataframename = outputfolder + "\\" + filename + "DataSheet" + "_ii" + str(integralint) + "_ct_" + str(correcthresh) + ".xlsx"
        microglialaveragename = outputfolder + "\\" + filename + "AveragesSheet" + "_ii" + str(integralint) + "_ct_" + str(correcthresh) + ".xlsx"
        # names

        MicroglialDataFrame.to_excel(microglialdataframename, index_label="Index")
        MicroglialAverages.to_excel(microglialaveragename, index_label="Index")
        # to excel

        print("cycle complete   |   " + str(round(((time.time() - currenttime) / 60), 3)) + " mins")

def imagedisplayer(inputfilelist, outputfolder, filepath, yeskey):

    #display images for editing (if they didn't choose to skip)

    will_restart = False
    firsttime = True
    counter = 0
    #set variables

    if yeskey[0]:

        listfiles = [image.path for image in os.scandir(outputfolder) if yeskey[1] in image.path]

    else:

        listfiles = [image.path for image in os.scandir(outputfolder)]

    #check whether a key was used. if so, slice for only those images with the key.

    if yeskey[2]:

        # if chose to skip editing

        for image, inputfile in zip(listfiles, inputfilelist):

            countera = 0
            name = ''

            while countera + 3 <= len(inputfile):

                if "TRI" == inputfile[countera:countera + 3]:

                    name = inputfile[countera + 3:]
                    break

                countera += 1

            shutil.copy(inputfile, "C:\pythonProject\HebertLabKM\Training\Imagesruncpij" + str(name))
            # copy files to imagej folder for sholl analysis

        print("\nskipped image edits...loading up ImageJ for sholl analysis!...\n")
        filemultiprocessing("C:\pythonProject\HebertLabKM\Training\Imagesruncpij")
        # load sholl analysis on that folder
        quit()

    else:

        # if not skipping editing

        for image, inputfile in zip(listfiles, inputfilelist):

            refimage = cv.imread(inputfile)
            prosimage = cv.imread(image)
            prosimage = cv.resize(prosimage, (int(screenWidth // 2) - 160, 1000))
            refimage = cv.resize(refimage, (int(screenWidth // 2) - 160, 1000))
            cv.imshow("prosimage", prosimage)
            cv.moveWindow('prosimage', 1419, 0)
            cv.imshow('ref', refimage)
            cv.moveWindow('ref', 0, 0)
            pyautogui.click(1980, 340)
            print("opened" + str(counter))
            # read/open images, and position them according to the center of the screen.

            while True:

                k = cv.waitKey(0) & 0xFF
                # indefintely pause until you get the right key input

                if k == 8:

                    if firsttime:

                        will_restart = True
                        # since an image has been edited, program will be run again on that image.
                        print("\nloading ImageJ...\n")
                        keycheck = input("\nDo you wish to add a keyword to the saved files? ('_copy' will be used as default) ('y' or 'n')...\n")

                        if keycheck in ['y', 'yes', 'yah', 'yeah', 'yup']:

                            keyword = input("\nenter the keyword...\n")

                        else:

                            keyword = "_copy"

                        firsttime = False
                        # store keyword to save copied images with ("copy" is used by default)

                    imagecopier(inputfilelist, counter, keyword)
                    cv.destroyAllWindows()
                    break

                # check for backspace ('no')--if entered, kill windows and load up imagej for editing. After editing is done, load next image.

                elif k == 13:

                    countera = 0
                    name = ''

                    while countera + 3 <= len(inputfile):

                        if "TRI" == inputfile[countera:countera + 3]:

                            name = inputfile[countera + 3:]
                            break

                        countera += 1

                    shutil.copy(inputfile, "C:\pythonProject\HebertLabKM\Training\Imagesruncpij" + str(name))
                    cv.destroyAllWindows()
                    break

                # check for the enter key ('ok'). if so, save image as a copy to the imagej folder for later, and go to the next image.

                elif k == 102:

                    flaggedimages.append(image)
                    print(image + " has been flagged\n")

                # check for f-key. means image has been flagged (will print a list of flagged images at the end).

                elif k == 103:

                    goodimages.append(image)
                    print(image + " has been marked as good\n")

                # check for the g-key. means image has been marked as good (will print a list of good images at the end).

            counter += 1
            # for indexing image files

        if will_restart:

            # will be true if you edited an image.
            print('\nprocessed all images. edits noted. re-running processing on modified images...\n')
            retain = input("\ndo you wish to retain settings?...\n")

            if retain in ['y', 'yes', 'yah', 'yeah', 'yup']:

                print("\nok, loading settings...\n")
                print('\nwarning: once you are done editing images and satisifed with cellprofiler output, quit, re-organize your input folder, and re-run all images...\n')
                main(True)

            else:

                print("\nsettings not retained\n")
                main(False)

            # note whether settings are retained, and re-run images. Be careful here--note what it's actually doing. It will simply re-run cellp on your new images (not all of the old ones). So, after you go through cycles of edits (and are happy with the cellp output on those edited images, compile all old good images + new images (deleting the ones you aren't using anymore), and re-run the entire program on them.
            # if you only run the program on an edited image, everything downstream will only run on that.

        else:

            print("\nno modifications made, loading up ImageJ for sholl analysis!...\n")
            print("\nhere were the flagged images from this experiment:\n")
            print(flaggedimages)
            print("\nhere were the good images from this experiment:\n")
            print(goodimages)
            filemultiprocessing("C:\pythonProject\HebertLabKM\Training\Imagesruncpij")
            # load up imagej sholl analysis
            quit()

def imagecopier(inputfilelist, counter, keyword):
    # make copies of images

    name = inputfilelist[counter]
    #get name based on index

    copyname = name[:len(name)-4] + keyword + ".tif"
    #create a copy name (need to remove tif off the end, then add keyword (to denote that it has been edited), then add it back on the end)

    shutil.copyfile(name, copyname)
    imageJeditor(copyname)
    #run opener on the copy

def imageJeditor(name):

    # open imageJ for editing images

    cv.destroyWindow('ref')
    ij = imagej.init(mode="interactive")
    ij.ui().showUI()
    # load imageJ and show the UI

    app = Application().connect(title_re = "ImageJ", found_index=0)
    dapp = app.top_window()
    dapp.move_window(screenWidth//4, int(screenHeight//2.5))
    # connect to the opened window (so we can move it around + make sure our clicks are accurate) + move image to the right position.

    pyautogui.click(689, 717)
    # click file
    pyautogui.click(724, 768)
    # click open
    time.sleep(1)
    pyautogui.write(name)
    # write the file name, pause so that first few characters aren't missed when writing.
    gx, gy = pyautogui.position()
    # get mouse position for reference
    pyautogui.click(gx + 325, gy + 29)
    # click to open file
    pyautogui.click(gx+ 795, gy - 100)
    # click on imageJ bar
    pyautogui.click(gx - 390, gy - 342)
    # click back on image
    dapp.move_window(0, screenHeight-400)
    # move imageJ bar below image
    pyautogui.click(217, 1317)
    # click on free draw
    print("\nhit 'enter' on your keyboard when you're done editing\n")
    keyboard.wait('enter')
    # wait for "enter" (i.e. that the file is good)
    pyautogui.click(58, 1279)
    # click file
    pyautogui.click(313, 1046)
    # click save
    pyautogui.click(1492, 534)
    # click to replace image in copied folder
    pyautogui.click(1357, 16)
    # click to close image
    ij.dispose()
    # close imageJ bar

def filemultiprocessing(folder):
    # multithreading so that tracking of elapsed time is possible (limit how long we try to do sholl analysis on a file/trace it), and so we don't hang
    # we need to track time elapsed here because imageJ has a tendency to hang (or take an unreasonable amount of time) on sholl analysis or tracing depending on image quality/how many cells there are.
    # what we need to do is save the traced files into a folder once it traces them, and then do automated sholl on them all at the end.
    # unfortunately, the only way to track whether they have been traced or not is to check a paramter via pyimageJ, which indicates whether the traces have been made. However, there is a delay from when these are made to when they are written to the file (and thus, when it can be saved as an actual "traces" file).
    # the delay is unpredictable, hence, I set it up in a way that it waits 60 seconds (can be edited, just will get less images the shorter you make it.) after which it quits the instance, and moves onto the next.
    # technically quitting the instance of imageJ unexpectedly would cause a normal program to error, but because of multiprocessing, we get to simply scrap that image (not included in sholl analysis), and start another instance/image.

    brokenlist = []
    # for storing which images ImageJ hangs on

    files = [file.path for file in os.scandir(folder)]
    counter = 0

    for file in files:

        if __name__ == '__main__':

            p = multiprocessing.Process(target=shollanalysis, args=(file, counter, files))
            p.start()

            if counter >= len(files)-1:

                while True:

                    try:

                        pd.read_csv("C:\pythonProject\HebertLabKM\Training\OutputData\_Sholl_Metrics.csv")
                        break

                    except:

                        pass

                time.sleep(900)
                p.terminate()

            # if we are on the last image, continously check whether sholl analysis has been completed. Once it has, wait 15 mins to give the file time to write.
            # you can defintely shorten the sleep time, but just give some delay for imagej/excel to write the file (otherwise it may be interrupted when writing, resulting in incomplete data).

            else:

                currenttime = time.time()

                while (time.time() - currenttime) < 300:

                    if not p.is_alive():

                        p.join()
                        break

            # if the program is dead (finished) before 300s, quit. If its been more than 300, move on to the next image.

                if p.is_alive():

                    print(file + " is broken")
                    brokenlist.append(file)
                    p.terminate()

                # if program is still alive after 300s (kill it, and mark the file as broken).

        counter += 1
        # for indexing purposes

    if len(brokenlist) > 0:

        print("the following files are bad: " + str(brokenlist))

    # if you have broken files, print them

    print("\ncalling data organizer...\n")
    dataorganizer()
    # call data organizer

def shollanalysis(file, counter, files):

    # open up imageJ and run the files through tracing and sholl analysis

    ij = imagej.init(["sc.fiji:fiji:2.14.0", "org.morphonets:SNT:4.2.1"], mode="interactive")
    # instantiate imagej

    SNTService = jimport('sc.fiji.snt.SNTService')
    SNTUI = jimport('sc.fiji.snt.SNTUI')

    sntService = ij.get(SNTService.class_)
    # Retrieve the SNTService instance

    imp = ij.py.to_imageplus(ij.io().open(file))
    imp.removeScale()
    # open imageJ ui and display image. remove scale (it affects the tracing/ sholl analysis otherwise).

    imp = ij.py.to_dataset(imp)
    # open file, convert to file and then a numpy readable format

    imp = ij.py.to_imageplus(imp[:, :, 1])
    ij.IJ.setAutoThreshold(imp, "Huang dark no-reset")
    # rgb image, get second channel (green) via slicing at 1.
    # threshold with auto Huang

    ij.IJ.run(imp, "Convert to Mask", "")  # apply the mask
    ij.IJ.saveAs(imp, "Tiff", "C:/pythonProject/HebertLabKM/Training/TRI/File_00d16c1f_2800_44ca_930a_decad9cf3ecd.tif (green).tif")
    # need to save mask as tiff, otherwise it won't work for tracing.

    curenttime = time.time()
    plugin = sntService.initialize(ij.py.to_imageplus(ij.io().open(file), True))
    # initialize snt

    time.sleep(.5)
    app = Application().connect(title_re="Run Auto-tracing?", found_index=0)
    dia = app["Run Auto-tracing?"]
    dia.move_window(0,506)
    time.sleep(.5)
    pyautogui.click(283, 760)
    # click auto-trace

    sapp = Application().connect(title_re="SNT v4.2.1", found_index=0)
    dia = sapp["SNT v4 2.1"]
    # connect to SNT tab

    while True:

        try:

            if sntService.getPaths():
                print(sntService.getPaths())
                time.sleep(80)
                break

        except:

            pass

    # check if there are paths--if there are, delay a little so that image is ready to save.

    dia.move_window(0, 0)
    time.sleep(1)
    # orient window
    pyautogui.click(52, 81)
    # click file
    time.sleep(1)
    pyautogui.click(254, 351)
    # click save as
    pyautogui.dragTo(618, 361, .5)
    # drag it to save, so menu is open
    pyautogui.click(618, 361)
    # click on save as
    time.sleep(.5)
    pyautogui.tripleClick(661, 975)
    # make sure you click-once doesn't highlight all text.
    pyautogui.press("backspace")
    # delete if something is there, just in case.
    time.sleep(.5)
    # delay to give time to write after deleting.
    pyautogui.write(file[:len(file) - 4] + ".traces")
    # write the file name, but cut off tiff and replace with traces.
    pyautogui.click(1143, 968)
    # click to save

    if counter >= len(files)-1:

        # if on the last image, run the sholl analysis

        currenttimet = time.time()
        time.sleep(2)
        print("nearly done")
        pyautogui.click(362, 75)
        time.sleep(1)
        pyautogui.click(362, 159)
        time.sleep(1)
        pyautogui.dragTo(987, 159)
        time.sleep(1)
        pyautogui.click(987,475)
        pyautogui.press('enter')
        # bullk sholl analysis (click happens really fast, so it seems like it's not clicking bulk sholl, but it is)

        while True:

            try:

                pd.read_csv("C:\pythonProject\HebertLabKM\Training\OutputData\_Sholl_Metrics.csv")
                time.sleep(600)
                break

            except:
                pass

        # wait till file path is present, and give time for writing.

    paths = sntService.getPaths()
    pyautogui.click(41, 73)
    time.sleep(1)
    pyautogui.click(93, 634)
    time.sleep(1)
    # quit

    try:

        app = Application().connect(title_re="Really Quit?", found_index=0)
        dia = app["Really Quit?"]
        dia.move_window(0, 506)
        pyautogui.click(176, 692)
        # click to quit
        time.sleep(.5)
        print("time: " + str(time.time() - curenttime) + "s paths: " + str(paths))

    except:

        pass

    # catch an error in quitting

    ij.dispose()
    # dispose imageJ

def dataorganizer(boolredo = False):
    # organize the data into a consistent, readable format (across both cellprofiler and imageJ)
    # can also do a double, cross-verification between cellprofiler and imageJ output, and save the results.

    shollsheet = pd.read_csv("C:\pythonProject\HebertLabKM\Training\OutputData\_Sholl_Metrics.csv")
    cellprofilersheet = pd.read_excel(outputfolder + "\\" + microglialdataframename)
    # pass in spreadsheet values/read files

    valuestokeep = {}
    shollcells = 0
    cellpcells = 0
    sholldictionary = {}
    cellpdictionary = {}
    # create cell counts and storage dictionaries to use in slicing

    indexvalues = []
    counter = 0

    for value in shollsheet['-']:

        indexvalues.append(counter)
        counter += 1

    shollsheet.insert(0, 'Index', indexvalues)
    # quick for-loop to count index values to append to the sholl sheet

    namelist = []
    for center, filename, endradius, index in zip(shollsheet['Center'], shollsheet['-'], shollsheet['End radius'], shollsheet['Index']):

        name = ''

        if not pd.isnull(endradius) and 15 < int(endradius) < 150:

            # filter for not having a null end radius (failed sholl), while also having a endradius between 15 and 150--to try to reduce how many false tracings are included

            try:

                count = 0

                while filename[count] != " ":

                    count += 1

                name = filename[:count]
                ccount = 0

                while center[ccount] != ",":

                    ccount += 1

                xvalue = int(center[:ccount])
                ccount += 1
                oldcount = ccount

                while center[ccount] != ",":

                    ccount += 1

                yvalue = int(center[oldcount + 1:ccount])

                if name in sholldictionary.keys():

                    sholldictionary[name].append([xvalue, yvalue, int(endradius), index])
                    valuestokeep[name].append([index, 0])

                else:

                    sholldictionary[name] = [[xvalue, yvalue, int(endradius), index]]
                    valuestokeep[name] = [[index, 0]]

                # store names/values

                shollcells += 1
                # increment cells

            except:
                print('cell thrown out')

        namelist.append(name)
        # append names

    # parse through sholl data and find the x and y values for the centroids.
    # assemble dictionary with filename: [x,y] structure.

    shollsheet['-'] = namelist
    # replace the overly-complex imageJ names with the normal file names.

    filenamelist = []
    cellprofilercells = 0
    for filename, xpos, ypos, mpos, index in zip(cellprofilersheet['FileName'], cellprofilersheet['MicrogliaXPosition'], cellprofilersheet['MicrogliaYPosition'], cellprofilersheet['MicrogliaMajorAxisLength'], cellprofilersheet['Index']):

        # basically doing the same organization for cellprofiler

        countr = 0
        countname = 0

        while countr < len(filename) - 4:

            if filename[countr:countr + 4] == "File":

                break

            countr += 1

        # get file names

        name = filename[countr:len(filename) - 4]
        filenamelist.append(name)
        if "." in name:

            print("this filename doesn't end  with 'tif': " + str(filename))
            sys.exit()

        # need to check for improperly constructed names

        try:

            if name in cellpdictionary.keys():

                cellpdictionary[name].append([xpos, ypos, mpos, index])

            else:

                cellpdictionary[name] = [[xpos, ypos, mpos, index]]

            cellprofilercells += 1

        except:

            print('falsevalue')

        # from the sheet, take the x and y values, along with major axis, and index.

    cellprofilersheet['FileName'] = filenamelist
    # replace filelist

    print("\npre-filter shollcellcount: " + str(shollcells))
    print("pre-filter cellpcellcount: " + str(cellprofilercells))
    # print out pre-filter, raw cell counts

    finalcellcounter = 0

    if booldefaultanalysis and not boolredo:

        # if using default values
        pixelrange = 1500
        axisrange = 1500

    # this is only relevant if using the cross-check (need to uncomment out the code below if you want to use)

    else:

        # if using custom values

        while True:

            try:

                pixelrange = float(input("\nwhat +/- for the centroid match range between cellp and sholl?...\n"))
                break

            except:

                print("\nyou entered an invalid value. please try again (enter a number)...\n")

        while True:

            try:

                axisrange = float(input("\nwhat +/- for the majoraxis(halved)/radius match range between cellp and sholl?...\n"))
                break

            except:

                print("\nyou entered an invalid value. please try again (enter a number)...\n")

        # take input for the value +/- ranges

    # this is only relevant if using the cross-check (need to uncomment out the code below if you want to use)

    ###for skey in sholldictionary.keys():
        ###for ckey in cellpdictionary.keys():
            ###if ckey == skey:
             ###   siterlist = sholldictionary[skey]
              ###  alreadymatched = []
               ### for val in siterlist:
                ###    sxval = val[0]
                 ###   syval = val[1]
                  ###  smval = val[2]
                   ### sindex = val[3]
                    ### citerlist = cellpdictionary[ckey]
                    ### for val in citerlist:
                     ###   if val in alreadymatched:
                      ###      pass
                       ### else:
                        ###    cxval = val[0]
                         ###   cyval = val[1]
                          ###  cmval = val[2]
                           ### cindex = val[3]
                            ### if cxval - pixelrange < sxval < cxval + pixelrange and cyval - pixelrange < syval < cyval + pixelrange and cmval / 2 - axisrange < smval < cmval / 2 + axisrange:
                             ###   if skey in valuestokeep.keys():
                              ###      valuestokeep[skey].append([sindex, cindex])
                               ###     finalcellcounter += 1
                                ###    alreadymatched.append(val)
                                  ###  break
                            ### else:
                            ###    valuestokeep[skey] = [[sindex, cindex]]
                            ###    finalcellcounter += 1
                            ###    alreadymatched.append(val)
                            ###    break
    # go through all cell x, and y. compare between sholl and cellp, if match, within range + major axis matches as well, add their indexes to a "values to keep" dictionary, which also stores the file names they are kept from.
    # run at different distances, and see where the values that acount for "already matched" and not diverge--that is when you should set it (so that its not doing a first come, first serve identification).

    sindexlist = []
    cindexlist = []
    # create index lists to append

    print("\nfinal image stats (applied to both cellp and sholl):\n")

    for key, val in zip(valuestokeep.keys(), valuestokeep.values()):

        print(str(key) + " has " + str(len(val)) + ' cells')

        for v in val:

            sindexlist.append(v[0])
            cindexlist.append(v[1])

    # append to separate keep lists and print some index-related data.

    shollsheet = shollsheet[shollsheet['Index'].isin(sindexlist)]
    shollsheet.rename(columns={'-': 'FileName'}, inplace=True)
    # slice the lists for only what matches the indexes.
    ### cellprofilersheet = cellprofilersheet[cellprofilersheet['Index'].isin(cindexlist)]
        # nothing to filter for cellprofiler if the cross-check is off

    print(str(len(valuestokeep.keys())) + " images")
    print('\ntotal cells in sholl: ' + str(len(sindexlist)))

    try:

        print("cells per image: " + str(shollcells / len(valuestokeep.values())))

    except:

        print("\nzero cells identified..redo with new params...\n")
        dataorganizer(True)

    # finally, show the counts.

    redobool = input("\nare you happy with the results? or would you like to redo with different params ('y'/'r')\n")

    if redobool in ['y', 'h', 'happy', 'yes']:

        print('ok, saving files...')
        shollsheet.to_excel("C:\pythonProject\HebertLabKM\Training\AnalysisOutput" + "\\" + "sholl_spliced" + str(microglialdataframename), index=False)
        cellprofilersheet.to_excel("C:\pythonProject\HebertLabKM\Training\AnalysisOutput" + "\\" + "cellp_spliced" + str(microglialdataframename), index=False)
        # save sliced sheets
        print('raw files saved, running downstream analysis...')
        downstream(cellprofilersheet, shollsheet)

    # save sliced sheets + run downstream analysis

    else:

        dataorganizer(True)

    # redo the analysis if unhappy.

def downstream(cellprofilersheet, shollsheet):

    # downstream analysis

    while True:

        try:

            msheet = input("\nEnter the mapping sheet\n")
            msheet = msheet[1:len(msheet) - 1]
            mappingsheet = pd.read_excel(msheet)
            break

        except:

            print("invalid path")

        # input mapping sheet for sorting files (note, all mapping sheets need to be formatted like the OG one Rohan gave me. I based it on that).

    CRETransplantFiles = []
    PLXTransplantFiles = []
    ControlFiles = []
    # create file lists for sorting

    mappingdictionary = {}
    for col1, col2 in zip(mappingsheet['OriginalName'], mappingsheet['AssignedName']):

        countr = 0

        while True:

            if col2[countr:countr + 5] == "File_":

                break

            countr += 1

        filename = col2[countr:len(col2) - 4]
        # slicing sheet for file names

        ### if "Tam" in col1 or "TAM" in col1:
            ### CRETransplantFiles.append(filename)
            # if you want to add this back, change the one below to "elif"

        if "GFP" in col1:

            PLXTransplantFiles.append(filename)

        elif "Iba1" in col1 or "iba1" in col1 or "IBA1" in col1 or "Tam" in col1 or "TAM" in col1:

            ControlFiles.append(filename)

        else:

            print(str(col1) + " failed ")
            print("your files are named incorrectly. check your key, quitting...")
            sys.exit()

        # check for keywords. replace this with whatever keywords you use.

    print(str(len(PLXTransplantFiles)) + " Number of Transplant files")
    print(str(len(ControlFiles)) + " Number of control files")

    shollsheetnames = ['End radius', 'Max inters. (sampled)', "Max inters. radius (sampled)",
                       'Sum inters. (sampled)', "Mean inters. (sampled)", "Ramification index (sampled)",
                       "Branching index (sampled)", "Enclosing radius (sampled)", "Intersecting radii",
                       "Sholl decay", "R^2", "r"]
    # these are the columns it'll calculate averages/run t-tests on.

    cellprofilernames = ['MicrogliaCD68MeanIntensity', 'MicrogliaCD68IntegratedIntensity', 'MicrogliaArea',
                         'MicrogliaMajorAxisLength', 'MajorAxis/MicrogliaArea', 'NucleiArea',
                         'NucleiCD68MeanIntensity', 'NucleiCD68IntegratedIntensity', 'NucleiArea/MicrogliaArea',
                         'NucleiXPosition', 'NucleiYPosition', 'MicrogliaXPosition', 'MicrogliaYPosition']
    # these are the columns it'll calculate averages/run t-tests on.

    CREsheet_sholl = shollsheet[shollsheet['FileName'].isin(CRETransplantFiles)]
    PLXsheet_sholl = shollsheet[shollsheet['FileName'].isin(PLXTransplantFiles)]
    Controlsheet_sholl = shollsheet[shollsheet['FileName'].isin(ControlFiles)]
    CREsheet_cellp = cellprofilersheet[cellprofilersheet['FileName'].isin(CRETransplantFiles)]
    PLXsheet_cellp = cellprofilersheet[cellprofilersheet['FileName'].isin(PLXTransplantFiles)]
    Controlsheet_cellp = cellprofilersheet[cellprofilersheet['FileName'].isin(ControlFiles)]
    # slice the sheets corresponding to the different file types.

    keylist = [pd.DataFrame.mean, pd.DataFrame.median, pd.DataFrame.std]
    colkeys = ['colmean', 'colmedian', 'colstd']
    treatgroups = ['Transplant', 'Control']
    # lists for functions/columns/treatmentgroups

    ShollSheetList = []
    CellPSheet = []
    # make lists for cellp sheet and sholl sheet

    ssheetcombined = {"TransplantSheetSholl": [PLXsheet_sholl, PLXTransplantFiles], "ControlSheetSholl": [Controlsheet_sholl, ControlFiles]}
    csheetcombined = {"TransplantSheetCellP": [PLXsheet_cellp, PLXTransplantFiles], "ControlSheetCellP": [Controlsheet_cellp, ControlFiles]}
    # combine sheets into dictionaries

    for sheet, sheetname in zip(ssheetcombined.values(), ssheetcombined.keys()):

        filemeanvalues = {'FileName': []}

        for filename in sheet[1]:

            tempsheet = sheet[0]

            if filename in list(tempsheet["FileName"]):

                filemeanvalues["FileName"].append(filename)
                tempsheet = tempsheet[tempsheet['FileName'] == filename]

                for column in tempsheet.columns:

                    if column in shollsheetnames:

                        meanvalue = pd.DataFrame.mean(tempsheet[column])

                        if column in filemeanvalues.keys():

                            filemeanvalues[column].append(meanvalue)

                        else:

                            filemeanvalues[column] = [meanvalue]

        ShollSheetList.append(pd.DataFrame(filemeanvalues))

    # calculate image-wise means for sholl analysis

    for sheet, sheetname in zip(csheetcombined.values(), csheetcombined.keys()):

        filemeanvalues = {'FileName': []}

        for filename in sheet[1]:

            tempsheet = sheet[0]

            if filename in list(tempsheet["FileName"]):

                filemeanvalues["FileName"].append(filename)
                tempsheet = tempsheet[tempsheet['FileName'] == filename]

                for column in tempsheet.columns:

                    if column in cellprofilernames:

                        meanvalue = pd.DataFrame.mean(tempsheet[column])

                        if column in filemeanvalues.keys():

                            filemeanvalues[column].append(meanvalue)

                        else:

                            filemeanvalues[column] = [meanvalue]

        CellPSheet.append(pd.DataFrame(filemeanvalues))

    # same calculations for cellprofiler

    CellPDict = {}
    ShollSDict = {}
    cgroupings = []
    sgroupings = []

    for val1, val2, name in zip(ShollSheetList, CellPSheet, treatgroups):

        ShollSDict[name] = val1
        CellPDict[name] = val2

    # create a dictionary to make references easier for the two sheets

    ShollDictionary = {}

    for function, keyname in zip(keylist, colkeys):

        for sheet in treatgroups:

            for name in shollsheetnames:

                column = ShollSDict[sheet][name]
                keyvalue = function(column)

                if name in ShollDictionary.keys():

                    ShollDictionary[name].append(keyvalue)

                else:

                    ShollDictionary[name] = [keyvalue]

            sgroupings.append(str(keyname) + str(sheet))

    # calculate across all images in a group (for sholl)

    CellPDictionary = {}

    for function, keyname in zip(keylist, colkeys):

        for sheet in treatgroups:

            for name in cellprofilernames:

                column = CellPDict[sheet][name]
                keyvalue = function(column)

                if name in CellPDictionary.keys():

                    CellPDictionary[name].append(keyvalue)

                else:

                    CellPDictionary[name] = [keyvalue]

            cgroupings.append(str(keyname) + str(sheet))

    # calculate across all images in a group (for cellp)

    Celldf = pd.DataFrame(CellPDictionary)
    Sholldf = pd.DataFrame(ShollDictionary)

    Celldf.insert(0, "TreatmentGroups", cgroupings)
    Sholldf.insert(0, "TreatmentGroups", sgroupings)
    # insert the treatment groups for easier readability

    Celldf.to_excel("C:/pythonProject/HebertLabKM/Training/AnalysisOutput" + "/" + "cellp_analyzed.xlsx", index=False)
    Sholldf.to_excel("C:/pythonProject/HebertLabKM/Training/AnalysisOutput" + "/" + "sholl_analyzed.xlsx", index=False)
    # to excel

    significant, nonsignificant = [], []
    cnonsignificantnames, snonsignificantnames = [], []
    # tracking sigifcance

    for sheet1, sheet2 in zip(treatgroups, treatgroups[1:]):

        if sheet1 != sheet2:

            for name in cellprofilernames:

                column1 = CellPDict[sheet1][name]
                column2 = CellPDict[sheet2][name]
                t_statistic, p_value = ttest_ind(column1, column2)

                if p_value <= 0.10:

                    significant.append("significant result: " + str(name) + " " + str(sheet1) + " against " + str(sheet2) + " | " + " p_value: " + str(p_value))

                else:

                    nonsignificant.append("not significant: " + str(name) + " " + str(sheet1) + " against " + str(sheet2) + " | " + " p_value: " + str(p_value))
                    cnonsignificantnames.append(name)

                # remove this logic "(if p>.10)" if you want to see all differences
                # currently will only print differences for comparisions with a p-value less than .10.
                # can see all differences in the sheet it prints out anyways.

    for sheet1, sheet2 in zip(treatgroups, treatgroups[1:]):

        if sheet1 != sheet2:

            for name in shollsheetnames:

                column1 = ShollSDict[sheet1][name]
                column2 = ShollSDict[sheet2][name]
                t_statistic, p_value = ttest_ind(column1, column2)

                if p_value <= 0.10:


                    significant.append("significant result: " + str(name) + " " + str(sheet1) + " against " + str(sheet2) + " | " + " p_value: " + str(p_value))

                else:

                    nonsignificant.append("not significant " + str(name) + " " + str(sheet1) + " against " + str(sheet2) + " | " + " p_value: " + str(p_value))
                    snonsignificantnames.append(name)

                # remove this logic "(if p>.10)" if you want to see all differences
                # currently will only print differences for comparisions with a p-value less than .10.
                # can see all differences in the sheet it prints out anyways.

    print("\n\nsignificant cellprofiler differences: \n")

    for column in [column for column in Celldf.columns[1:] if column not in cnonsignificantnames]:

        index = 0

        for value in Celldf[column]:

            if index == 0:

                subvalue = value

            elif index == 1:

                if float(subvalue) - float(value) > 0:

                    operator = "+"

                else:

                    operator = ""

                print(str(column) + " difference for PLX-Control: " + operator + str(subvalue - value))
                break

            index += 1

    # print significant cellp differences

    print("\n\nsignificant sholl differences: \n")

    for column in [column for column in Sholldf.columns[1:] if column not in snonsignificantnames]:

        index = 0

        for value in Sholldf[column]:

            if index == 0:

                subvalue = value

            elif index == 1:

                if float(subvalue) - float(value) > 0:

                    operator = "+"

                else:

                    operator = ""

                print(str(column) + " difference for PLX-Control: " + operator + str(subvalue - value))

            index += 1

    # print significant sholl differences

    print("\n\n\nt-test results...\n")
    print("Significant Results...\n")

    for value in significant:

        print(value)

    print("\nNon-significant Results...\n")

    for value in nonsignificant:

        print(value)

    # print p-values

def imageJtcoordtracer(file):
    # coordinate tracing for density analysis

    counter = len(file)

    while file[counter-1:counter] != "\\":

        counter-=1

    # get image name

    filename = file[counter:len(file)]
    ij = imagej.init(mode="interactive")
    ij.ui().showUI()
    #read/open imageJ
    app = Application().connect(title_re="ImageJ", found_index=0)
    dapp = app.top_window()
    #connect to the top window of imageJ
    imp = ij.IJ.openImage(file)
    imp.removeScale()
    #open file and remove scale. KEY: must open with "openImage", otherwise, .io().open() will return the wrong image type.
    ij.ui().show(imp)
    #show image
    dapp.move_window(0, screenHeight - 400)
    time.sleep(.1)
    #need to add sleep delays so it can find the window (if really wanted to, could do try/except loop)
    app = Application().connect(title_re=filename, found_index=0)
    dapp = app.top_window()
    #connect to image file
    imagewidth = dapp.get_properties()['rectangle'].right - dapp.get_properties()['rectangle'].left
    dapp.move_window(screenWidth//2-imagewidth//2, 0)
    time.sleep(.5)
    #find the width of the image, and use that to move it into right place (so your click works)
    pyautogui.rightClick(270,1319)
    pyautogui.click(277, 1368)
    pyautogui.click(screenWidth//2, 268)
    #click on free draw, and then click back on image

    print("\ntracing order: 1. Top of Neocortex 2. Bottom of Neocortex 3. Top of Hippocampus 4. Bottom of Hippocampus\n")
    print("after tracing a section, press 'b', and don't touch the mouse, until the clicks stop\n")
    counter = 0

    while counter < 4:

        keyword = "NeocortexPoints"

        if counter >= 2:

            keyword = "HippocampusPoints"

        if counter % 2 == 0:

            orientation = "Top"

        else:

            orientation = "Bottom"

        # track what part you are tracing

        name = filename[:len(filename)-4] + keyword + orientation
        keyboard.wait("b")
        # wait for it to b input (will open up the points), and then send commands to keyboard.
        keyboard.press("y")
        # y is a hotkey in imageJ
        time.sleep(.8)
        app = Application().connect(title_re="Properties", found_index=0)
        dapp = app.top_window()
        #connect to properties tab
        prowidth = dapp.get_properties()['rectangle'].right - dapp.get_properties()['rectangle'].left
        dapp.move_window(screenWidth // 2 - prowidth // 2, 780)
        # get the window properties, and move to the right location
        pyautogui.click(990, 1193)
        pyautogui.press("enter")
        #click to get coordinates, and enter to proceed
        time.sleep(.1)
        pyautogui.hotkey('ctrl', 's')
        #click to save file
        time.sleep(.5)
        #this delay is crucial for keyboard writing time.
        pyautogui.write("C:\\pythonProject\\HebertLabKM\\Training\\DensityCoords\\" + name + ".csv")
        pyautogui.press("enter")
        #save filename
        time.sleep(.1)
        app = Application().connect(title_re=filename, found_index=0)
        dapp = app.top_window()
        dapp.close()
        #find the coordinates window, and close it
        counter+=1

    ij.dispose()
    #dispose imageJ

def densitypreprocessing(cellprofilerfile):
    # analagous to data organizer, but for density counts

    cellprofilersheet = pd.read_excel(cellprofilerfile)
    cellprofilercells = 0
    filenamelist = []
    uniquenamelist = []
    # using this, because set function doesn't guarantee order preservation

    for filename, xpos, ypos in zip(cellprofilersheet['FileName'], cellprofilersheet['MicrogliaXPosition'], cellprofilersheet['MicrogliaYPosition'],):

        countr = len(filename)

        while filename[countr - 1:countr] != "\\":

            countr -= 1

        name = filename[countr:len(filename) - 4]
        filenamelist.append(name)
        cellprofilercells += 1

        if name not in uniquenamelist:

            uniquenamelist.append(name)

    # get the right file names

    cellprofilersheet['FileName'] = filenamelist
    # replace names

    print("avg pre-filter summed microglialcellcount: " + str(int(cellprofilercells/len(uniquenamelist))) + "\n")
    # print out pre-filter, raw cell counts

    cellprofilersheet.to_excel("C:/pythonProject/HebertLabKM/Training/AnalysisOutput" + "/" + "cellp_spliced" + str(microglialdataframename), index=False)
    # save file

    for filename in uniquenamelist:

        print("\n\n" + filename)
        densitycounter(filename, cellprofilersheet)

    # iterate through names, and call density counter on them

def densitycounter(reffilename, cellprofilersheet):

    densityfiles = [file.path for file in os.scandir("C:\pythonProject\HebertLabKM\Training\DensityCoords") if reffilename in file.path]
    NeocortexFiles = [file for file in densityfiles if "Neocortex" in file]
    HippocampusFiles = [file for file in densityfiles if "Hippocampus" in file]
    # list comprehension for files

    HippocampusDictionary, NeocortexDictionary = {}, {}
    mediansum = 0
    counter = 0
    hippomax, neomax = 0, 0
    filelist = [NeocortexFiles, HippocampusFiles]
    # variables

    for refname in filelist:

        for fname in refname:

            referencename = fname

            if "Top" in fname:

                refnametop = pd.read_csv(fname)

            elif "Bottom" in fname:

                refnamebottom = pd.read_csv(fname)

            # read the corresponding coordinate file based on what you're processing

        startmedx = (pd.DataFrame.median(refnametop['X']) + pd.DataFrame.median(refnamebottom['X']))//2
        mediansum += startmedx
        # for calculating the image median

        max_x = max(pd.DataFrame.max(refnametop['X']), pd.DataFrame.max(refnamebottom['X']))
        min_x = min(pd.DataFrame.min(refnametop['X']), pd.DataFrame.min(refnamebottom['X']))
        # max and minimum x's

        increment = 10
        # how many pixels wide the integration is

        if "Hippocampus" in referencename:

            hippomax = max_x

        else:

            neomax = max_x

        # store max values

        curpos = max_x
        # start at the maximum x value (work down from there, looking for cells)

        while curpos >= min_x + increment:

            for filename, xpos, ypos, index in zip(cellprofilersheet['FileName'], cellprofilersheet['MicrogliaXPosition'], cellprofilersheet['MicrogliaYPosition'], cellprofilersheet['Index']):

                if filename in reffilename and curpos - increment < xpos <= curpos:
                    # if cell is within increment at current tposition

                    firstxmintop = [100000, -10000]
                    closestdistancetop = 0

                    for xvalue, yvalue in zip(refnametop['X'], refnametop['Y']):

                        if abs(xvalue - xpos) <= abs(firstxmintop[0] - xpos):

                            if abs(xvalue - xpos) == abs(firstxmintop[0] - xpos):

                                if yvalue > firstxmintop[1]:

                                    firstxmintop = [xvalue, yvalue]
                                    closestdistancetop = abs(firstxmintop[0] - xpos)
                            else:

                                firstxmintop = [xvalue, yvalue]
                                closestdistancetop = abs(firstxmintop[0] - xpos)

                    # this is to calculate top closest x values and y values of the coordinate tracing to check whether cell fits between them

                    firstxminbottom = [100000, 100000]
                    closestdistancebottom = 0
                    for xvalue, yvalue in zip(refnamebottom['X'], refnamebottom['Y']):

                        if abs(xvalue - xpos) <= abs(firstxminbottom[0] - xpos):

                            if abs(xvalue - xpos) == abs(firstxminbottom[0] - xpos):

                                if yvalue < firstxminbottom[1]:

                                    firstxminbottom = [xvalue, yvalue]
                                    closestdistancebottom = abs(firstxminbottom[0] - xpos)

                            else:

                                firstxminbottom = [xvalue, yvalue]
                                closestdistancebottom = abs(firstxminbottom[0] - xpos)

                    # this is to calculate bottom closest x values and y values of the coordinate tracing to check whether cell fits between them

                    closeyvalues = [firstxmintop[1], firstxminbottom[1]]
                    closestdistance = (closestdistancetop + closestdistancebottom)//2
                    # get the closest y values, and average closest distance value (across bottom and top reference points)

                    if closestdistance <= 100 and min(closeyvalues) <= ypos <= max(closeyvalues):
                        # make sure closest values are actually close (so they are relevant as boundaries) + check that the cell fits between the minimum and maximum.

                        area = (max(closeyvalues) - min(closeyvalues)) * (increment)
                        # area for density (approximates a rectangle)

                        if "Hippocampus" in referencename:

                            if filename in HippocampusDictionary.keys():

                                HippocampusDictionary[filename].append([index, curpos - increment, area, curpos])

                            else:

                                HippocampusDictionary[filename] = [[index, curpos - increment, area, curpos]]

                        else:

                            if filename in NeocortexDictionary.keys():

                                NeocortexDictionary[filename].append([index, curpos - increment, area, curpos])

                            else:

                                NeocortexDictionary[filename] = [[index, curpos - increment, area, curpos]]

                        # append index, minpos, area, and maxposition to appropriate dictionary

            curpos -= increment
            # increment

    avgmed = mediansum//2
    print("\napprox median of image: " + str(avgmed))
    # average median (between neocortex and hippocampal files)

    for key in NeocortexDictionary.keys():

        print("\nneocortical cell count: " + str(len(NeocortexDictionary[key])))
        curposlist = [neomax]

        for keyval in NeocortexDictionary[key]:

            curposlist.append(keyval[3])
            # store positions

        curposlist = sorted(set(curposlist))
        # sort the list and use set function
        greatestdiff = 0
        avggap = 0
        lastvalue = min(curposlist)
        # end value

        for value in curposlist:

            if value-lastvalue > greatestdiff:

                greatestdiff = value-lastvalue
                rightval = value
                leftval = lastvalue

            avggap += value-lastvalue
            lastvalue = value

        avggap = avggap / (len(curposlist)-1)
        print("greatest Neocortical gap: " + str(greatestdiff) + "; starting at " + str(leftval) + ", until " + str(rightval))
        print("avg Neocortical gap: " + str(avggap))
        # greatest gap/average neocortical gap calculated

    for key in HippocampusDictionary.keys():

        print("\nhippocampal cell count: " + str(len(HippocampusDictionary[key])))
        curposlist = [hippomax]

        for keyval in HippocampusDictionary[key]:

            curposlist.append(keyval[3])
        # store positions

        curposlist = sorted(set(curposlist))
        # sort the list and use set function
        greatestdiff = 0
        avggap = 0
        lastvalue = min(curposlist)
        # end value

        for value in curposlist:

            if value - lastvalue > greatestdiff:

                greatestdiff = value - lastvalue
                rightval = value

                leftval = lastvalue

            avggap += value - lastvalue
            lastvalue = value

        avggap = avggap / (len(curposlist)-1)
        print("greatest Hippocampal gap: " + str(greatestdiff) + "; starting at " + str(leftval) + ", until " + str(rightval))
        print("avg Hippocampal gap: " + str(avggap))
        # greatest gap/average neocortical gap calculated

    curposdictionary = {'x-range': [], 'cell_count': [], 'count_density': [], 'sample_area': []}
    for key in HippocampusDictionary.keys():
        # compress everything into a dataframe/also, if you can get good measurements from the low zoom images (I doubt it), I included functionality for calculating the means across the cells at each position range.

        valuelist = []
        throwaway = [valuelist.append([value[1], value[2], value[3]]) for value in HippocampusDictionary[key] if [value[1], value[2], value[3]] not in valuelist]
        # use a comprehension for appending to get linear evaluation time

        curpos = min([value[1] for value in HippocampusDictionary[key]])
        while curpos <= max_x:

            for value in valuelist:

                if value[0] == curpos:

                    area = value[1]
                    finalpos = value[2]
                    break

            indexvalues = [value[0] for value in HippocampusDictionary[key] if value[1] == curpos]
            tempsubset = cellprofilersheet[cellprofilersheet['Index'].isin(indexvalues)]

            if len(indexvalues) > 0:

                for column in tempsubset.columns:

                    if column not in ["FileName", "Index", "MicrogliaXPosition", "MicrogliaYPosition"]:

                        mean = pd.DataFrame.mean(tempsubset[column])

                        if column in curposdictionary.keys():

                            curposdictionary[column].append(mean)

                        else:

                            curposdictionary[column] = [mean]

                curposdictionary['cell_count'].append(len(indexvalues))
                curposdictionary['count_density'].append(len(indexvalues) / area)
                curposdictionary['x-range'].append(str(int(min(curpos, finalpos))) + ' : ' + str(int(max(curpos, finalpos))))
                curposdictionary['sample_area'].append(area)

            curpos += increment

        pd.DataFrame(curposdictionary).to_excel("C:/pythonProject/HebertLabKM/Training/DensityOutput" + "/" + str(key) + "_" + str(avgmed) + "_" + "HippocampalDistributiton" + ".xlsx", index=False)

    curposdictionary = {'x-range': [], 'cell_count': [], 'count_density': [], 'sample_area': []}
    for key in NeocortexDictionary.keys():
        # compress everything into a dataframe/also, if you can get good measurements from the low zoom images (I doubt it), I included functionality for calculating the means across the cells at each position range.

        valuelist = []
        throwaway = [valuelist.append([value[1], value[2], value[3]]) for value in NeocortexDictionary[key] if [value[1], value[2], value[3]] not in valuelist]
        # use a comprehension for appending to get linear evaluation time

        curpos = min([value[1] for value in NeocortexDictionary[key]])
        while curpos <= max_x:

            for value in valuelist:

                if value[0] == curpos:

                    area = value[1]
                    finalpos = value[2]
                    break

            indexvalues = [value[0] for value in NeocortexDictionary[key] if value[1] == curpos]
            tempsubset = cellprofilersheet[cellprofilersheet['Index'].isin(indexvalues)]

            if len(indexvalues) > 0:

                for column in tempsubset.columns:

                    if column not in ["FileName", "Index", "MicrogliaXPosition", "MicrogliaYPosition"]:

                        mean = pd.DataFrame.mean(tempsubset[column])

                        if column in curposdictionary.keys():

                            curposdictionary[column].append(mean)

                        else:

                            curposdictionary[column] = [mean]

                curposdictionary['cell_count'].append(len(indexvalues))
                curposdictionary['count_density'].append(len(indexvalues) / area)
                curposdictionary['x-range'].append(str(int(min(curpos, finalpos))) + ' : ' + str(int(max(curpos, finalpos))))
                curposdictionary['sample_area'].append(area)

            curpos += increment

        pd.DataFrame(curposdictionary).to_excel("C:/pythonProject/HebertLabKM/Training/DensityOutput" + "/" + str(key) + "_" + str(avgmed) + "_" + "NeocorticalDistribution" + ".xlsx", index=False)

def main(x=False):
    # this is the main function (hence its name). It is the only function that is directly called (all other functions are called from calls to other functions)
    # unless you absolutely need to, avoid editing this function. It utilizes multiprocessing so that it can create instances of the jvm for cellprofiler and imageJ, without getting an error for using the jvm multipurposely.
    # if you do need to change it, read up on multiprocessing to avoid unintended side effects. The evaluation order is confusing.

    if __name__ == "__main__":

        # you cannot call main from another file

        filelist = inputfunction(x)
        # call the input function, and store it as "filelist"

        p = Process(target=cellprofiler, args=(filelist[0], filelist[2], filelist[3], filelist[4], filelist[5], filelist[6], filelist[8]))
        p.start()
        p.join()
        # run cellprofiler

        if not filelist[8]:

            imagedisplayer(filelist[0], filelist[3], filelist[1], filelist[7])

        else:

            for image in filelist[0]:

                imageJtcoordtracer(image)

            densitypreprocessing(outputfolder + "\\" + microglialdataframename)

        # run functions downstream

        print("\nanalysis is complete")

main()