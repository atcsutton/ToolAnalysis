########################################################################################################################
#
# 8888888b.  d8b                         .d8888b.                             888    d8b
# 888   Y88b Y8P                        d88P  Y88b                            888    Y8P
# 888    888                            888    888                            888
# 888   d88P 888 88888b.   .d88b.       888         .d88b.  888  888 88888b.  888888 888 88888b.   .d88b.
# 8888888P"  888 888 "88b d88P"88b      888        d88""88b 888  888 888 "88b 888    888 888 "88b d88P"88b
# 888 T88b   888 888  888 888  888      888    888 888  888 888  888 888  888 888    888 888  888 888  888
# 888  T88b  888 888  888 Y88b 888      Y88b  d88P Y88..88P Y88b 888 888  888 Y88b.  888 888  888 Y88b 888
# 888   T88b 888 888  888  "Y88888       "Y8888P"   "Y88P"   "Y88888 888  888  "Y888 888 888  888  "Y88888
#                              888                                                                     888
#                         Y8b d88P                                                                Y8b d88P
#                          "Y88P"                                                                  "Y88P"
# --------------------------------------------------------------------------------------------------------------------
#                  Daniel Tobias Schmid, Feb. 2023, dschmid@fnal.gov / d.schmid@students.uni-mainz.de
# --------------------------------------------------------------------------------------------------------------------
#
# The Cherenkov-ring-counting tool is used to classify events as single- or multi-ring by analyzing 10x16 PMT hit maps.
#   Data is loaded from a CSV file of CNNImage displays, or using the CNNImage tool can be created with data provided
#   by other tools. It uses a keras/tf based Convolutional Neural Network. To use this tool in a ToolChain, users
#   must populate the the corresponding "/configfiles/YourToolChain/RingCountingConfig" file with the following
#   information:
#     (Square brackets are a variable in the config file)
#     1. The path to a file containing a list of files of the PMT data to be used, if data is to be loaded from a csv
#        file. Must be in CNNImage format.
#        -> defined by setting [[files_to_load]]
#     2. Whether to load from a csv file or use a previous ToolChain to process/load processed data/MC
#        -> defined by setting [[load_from_csv]]
#     3. Whether to save to a csv file when previously loading data from a csv file (both have to be 1)
#        -> defined by setting [[save_to_csv]]
#     4. The model version
#        -> defined by setting [[version]]
#     5. The model directory path (not including the model's filename, but including last slash)
#        -> defined by setting [[model_path]]
#     6. Which PMT mask to use (some PMTs have been turned off in the training); check documentation for which model
#        requires what mask.
#       -> defined by setting [[pmt_mask]]
#     7. Where to save the predictions, when save_to_csv == 1
#       -> defined by setting [[save_to]]
#     8. Whether a single model or ensemble should be used
#       -> defined by setting [[model_is_ensemble]]
#     9. How many models make up the ensemble. If model count is N, "sub"-models are labeled 0, 1, ..., N-1.
#       -> defined by setting [[ensemble_model_count]]
#     10. How the model predictions should be combined when using an ensemble. Supported:
#         - "None" (the type), only the first model's predictions are used. (blank line in config file)
#         - "average", average predictions of all models
#         - "voting", average predictions and in addition a majority-voting prediction is also produced.
#       -> defined by setting [[ensemble_prediction_combination_mode]]
#
#   An example config file can be found in the configfiles/RingCounting/ ToolChain.
#
#
#  When using on the grid, make sure to only use onsite computing resources. TensorFlow is not supported at all offsite
#    facilities (July 2024). Also make sure to send the model to the grid, bundled with your process.
#
#
# Documentation of the tool, model versions and performance can(will) be found at (anniegpvm-machine):
#     /pnfs/annie/persistent/users/dschmid/RingCountingStore/documentation/ **TODO**
# All models are located in (anniegpvm-machine):
#     /pnfs/annie/persistent/users/dschmid/RingCountingStore/models/
#
########################################################################################################################
from Tool import *

import numpy as np
import tensorflow as tf


class RingCountingGlobals:
    # The reason to mask some PMTs to 0 is to disable PMTs that differ in the MC and experimental datasets.
    #   If a model has been trained with a specific mask toggled on, it is crucial that when evaluating experimental
    #   data with this model, that the same specific mask is also toggled on. This is due to the model not having
    #   encountered these PMTs being enabled during training, and thus probably producing worse predictions with
    #   the un-encountered PMTs enabled when predicting.
    PMT_MASKS = {
        "none": [],
        # From per-PMT p.e. distribution histograms -> which PMTs have different MC / data curves
        "pe_curve_divergent_november_22": [56, 58, 84, 85, 87, 88, 89, 100, 118, 120, 121, 137, 138, 139],
        # Difference in MC and data on specific datasets used for training
        "active_divergent_november_22": [22, 124],
        # Both november 22 lists
        "november_22": [22, 56, 58, 84, 85, 87, 88, 89, 100, 118, 120, 121, 124, 137, 138, 139],
        # Load from config TODO
        "custom": []
    }


class RingCounting(Tool, RingCountingGlobals):
    # ----------------------------------------------------------------------------------------------------
    # Data stuff
    cnn_image_pmt = None  # To be loaded with PMT data

    # ----------------------------------------------------------------------------------------------------
    # Config stuff
    load_from_csv = std.string()  # if 1, load 1 or more CNNImage formatted csv file instead of using toolchain
    save_to_csv = std.string()  # if 1, save as a csv file in format MR prediction, SR prediction
    files_to_load = std.string()  # List of files to be loaded (must be in CNNImage format,
    #   load_from_csv has to be true)
    version = std.string()  # Model version
    model_path = std.string()  # Path to model directory
    pmt_mask = std.string()  # See RingCountingGlobals
    save_to = std.string()  # Where to save the predictions to
    model_is_ensemble = std.string()  # Whether the model consists of multiple models acting as a mixture of experts
    #   (MOE)/ensemble
    ensemble_model_count = std.string()  # Count of models used in the ensemble
    ensemble_prediction_combination_mode = std.string()  # How predictions of models are combined: average, voting, ..

    # ----------------------------------------------------------------------------------------------------
    # Model variables
    model = None  # Union[TF.model/Keras.model, None]
    ensemble_models = None  # Union[List[TF.model/Keras.model], None]
    predicted = None  # np.array()
    predicted_ensemble = None  # List[np.array()]

    def Initialise(self):
        """ Initialise RingCounting tool object in following these steps:
        1. Load necessary config
        2. Load data
          2.1 Modify data
        3. Load model

        For specific details, see compartmentalised function docstrings.
        """
        # ----------------------------------------------------------------------------------------------------
        # Debug area
        self.m_log.Log(__file__ + " Initialising", self.v_debug, self.m_verbosity)

        # ----------------------------------------------------------------------------------------------------
        # Config area
        self.m_variables.Get("files_to_load", self.files_to_load)
        self.files_to_load = str(self.files_to_load)  # cast to str since std.string =/= str
        self.m_variables.Get("load_from_csv", self.load_from_csv)
        self.load_from_csv = "1" == str(self.load_from_csv)
        self.m_variables.Get("save_to_csv", self.save_to_csv)
        self.save_to_csv = "1" == str(self.save_to_csv)
        self.m_variables.Get("version", self.version)
        self.m_variables.Get("model_path", self.model_path)
        self.m_variables.Get("pmt_mask", self.pmt_mask)
        self.m_variables.Get("save_to", self.save_to)
        self.save_to = str(self.save_to)  # cast to str since std.string =/= str
        self.pmt_mask = self.PMT_MASKS[self.pmt_mask]
        self.m_variables.Get("model_is_ensemble", self.model_is_ensemble)
        self.model_is_ensemble = "1" == self.model_is_ensemble
        self.m_variables.Get("ensemble_model_count", self.ensemble_model_count)
        self.ensemble_model_count = int(self.ensemble_model_count)
        if self.ensemble_model_count % 2 == 0:
            self.m_log.Log(__file__ + f" WARNING: Number of models in ensemble is even"
                                      f" ({self.ensemble_model_count}). Can lead to unexpected classification when"
                                      f" using voting to determine ensemble predictions.",
                           self.v_warning, self.m_verbosity)
        self.m_variables.Get("ensemble_prediction_combination_mode", self.ensemble_prediction_combination_mode)
        if self.ensemble_prediction_combination_mode not in [None, "average", "voting"]:
            self.m_log.Log(__file__ + f" WARNING: Unsupported prediction combination mode selected"
                                      f" ({self.ensemble_prediction_combination_mode}). Defaulting to 'average'.",
                           self.v_warning, self.m_verbosity)
            self.ensemble_prediction_combination_mode = "average"

        # ----------------------------------------------------------------------------------------------------
        # Loading data
        if self.load_from_csv:
            self.load_data()
        else:
            self.m_log.Log(__file__ + " Not loading data from csv file.", self.v_message, self.m_verbosity)
            self.cnn_image_pmt = np.array([])

        # ----------------------------------------------------------------------------------------------------
        # Loading model
        self.load_model()

        return 1

    def Execute(self):
        """ Execute the tool by generating model predictions on the supplied data. """
        self.m_log.Log(__file__ + " Executing", self.v_debug, self.m_verbosity)

        self.get_next_event()
        self.mask_pmts()
        self.predict()

        self.process_predictions()

        return 1

    def Finalise(self):
        """ Finalise the tool by saving the predictions. """
        self.m_log.Log(__file__ + " Finalising", self.v_debug, self.m_verbosity)
        if self.save_to_csv and self.load_from_csv:
            # TODO: ToolChain -> save to csv
            # Can currently only do csv -> csv or ToolChain -> BoostStore.
            self.save_data()

        return 1

    def load_data(self):
        """ Load data in the CNNImage format.

        The files to load are specified in a txt file which filepath is defined in the RingCountingConfig file under
        [[files_to_load]]. The CNNImage format is defined in the CNNImage tool. In short, a 160-entry list of numbers
        (cast-able to numpy float-types) in comma separated format, where events are separated by newlines, is expected.

        Uncommenting a filepath in the [[files_to_load]] file is supported by placing "#" as the first character.

        Handles empty files by skipping and logging a warning.
        """
        self.cnn_image_pmt = np.array([])
        # Removes trailing and leading whitespace and "\n"
        self.files_to_load = self.files_to_load.strip(" \n")
        with open(self.files_to_load, "r") as file:
            for data_path in file:
                # Removes trailing and leading whitespace and "\n"
                data_path = data_path.strip(" \n")
                # Supports commenting out a file
                if data_path[0] == "#":
                    continue

                arr = np.loadtxt(data_path, delimiter=",")
                if arr.shape == (160,):
                    # Explicitly adding the extra dimension using np.array([]) is done to ensure the data is
                    #   reshaped into the shape (-1, 160) in case of having only a single entry.
                    arr = np.array([arr])
                arr = np.reshape(arr, (-1, 160))

                if len(self.cnn_image_pmt) > 0:
                    self.cnn_image_pmt = np.concatenate([arr, self.cnn_image_pmt], axis=0)
                    self.m_log.Log(__file__ + f" Successfully loaded CNNImage file {data_path} with {len(arr)} events.",
                                   self.v_debug, self.m_verbosity)
                elif len(arr) > 0:
                    self.m_log.Log(__file__ + f" Successfully loaded CNNImage file {data_path} with {len(arr)} events.",
                                   self.v_debug, self.m_verbosity)
                    self.cnn_image_pmt = arr
                else:
                    self.m_log.Log("WARNING: " + __file__ + f" Attempted to load an empty PMT-datafile (CNNImage). "
                                                            f"({data_path})",
                                   self.v_debug, self.m_verbosity)

    def save_data(self):
        """ Save the data to the specified [[save_to]]-file. When using an ensemble, each line contains all of the
        individual model's predictions for that event (ordered as MR1,SR1,MR2,SR2,...).
        """
        if self.model_is_ensemble:
            np.savetxt(self.save_to, np.array(self.predicted_ensemble).flatten(), delimiter=",")
        else:
            np.savetxt(self.save_to, self.predicted, delimiter=",")

    def mask_pmts(self):
        """ Mask PMTs to 0. The PMTs to be masked is given as a list of indices, defined by setting [[pmt_mask]].
        For further details check the RingCountingGlobals class.
        """
        if self.load_from_csv:
            for event in self.cnn_image_pmt:
                np.put(event, self.pmt_mask, 0, mode='raise')
        else:
            np.put(self.cnn_image_pmt, self.pmt_mask, 0, mode='raise')

    def load_model(self):
        """ Load the specified model [[version]]. If [[model_is_ensemble]], load all models in ensemble.
        Models files are expected to be named as 'model_path + RC_model_v[[version]].model' for single models, and
        'model_path + RC_model_ENS_v[[version]].i.model', where i in {0, 1, ..., [[ensemble_model_count]] - 1} for
        ensemble models.
        """
        if self.model_is_ensemble:
            self.ensemble_models = [
                tf.keras.models.load_model(self.model_path + f"RC_model_ENS_v{self.version}.{i}.model")
                for i in range(0, self.ensemble_model_count)
            ]
        else:
            self.model = tf.keras.models.load_model(self.model_path + f"RC_model_v{self.version}.model")

    def get_next_event(self):
        """ Get the next event from the BoostStore. """
        if self.load_from_csv:
            return

        reco_event_bs = self.m_data.Stores.at("RecoEvent")
        get_ok = reco_event_bs.Has("CNNImageCharge")
        self.cnn_image_pmt = std.vector['double'](range(160))

        if get_ok:
            reco_event_bs.Get("CNNImageCharge", self.cnn_image_pmt)
        else:
            self.m_log.Log(__file__ + " ERROR: CNNImageCharge not present in RecoEvent boost store.",
                           self.v_error, self.m_verbosity)

        # loop over std::vector to convert to list
        self.cnn_image_pmt = [x for x in self.cnn_image_pmt]
        # Explicitly adding the extra dimension using np.array([]) is done to ensure the data is properly
        #   reshaped into the shape (-1, 160).
        self.cnn_image_pmt = np.array([self.cnn_image_pmt])
        self.cnn_image_pmt = np.reshape(self.cnn_image_pmt, (-1, 160))

    def predict(self):
        """
        Classify events in single- and multi-ring events using a keras model. Store a list of 2-dimensional predictions
        (same order as input) to self.predicted. Predictions are given as [MR prediction, SR prediction].
        """

        self.m_log.Log(__file__ + " PREDICTING", self.v_message, self.m_verbosity)
        if self.model_is_ensemble:
            self.predicted_ensemble = [
                m.predict(np.reshape(self.cnn_image_pmt, newshape=(-1, 10, 16, 1))) for m in self.ensemble_models
            ]
        else:
            self.predicted = self.model.predict(np.reshape(self.cnn_image_pmt, newshape=(-1, 10, 16, 1)))

    def process_predictions(self):
        """ Process the model predictions. If an ensemble is used, calculate final predictions based on the selected
        ensemble mode. Finally, store predictions in the RecoEvent BoostStore.

        Store the output of the averaging ensemble and single model within the RecoEvent BoostStore under
          RingCountingSRPrediction,
          RingCountingMRPrediction.
        Store the voted-for class prediction of the voting ensemble within the RecoEvent BoostStore under
          RingCountingVotingSRPrediction,
          RingCountingVotingMRPrediction.
        """
        predicted_sr = -1
        predicted_mr = -1
        reco_event_bs = self.m_data.Stores.at("RecoEvent")

        if self.model_is_ensemble:
            if self.ensemble_prediction_combination_mode is None:
                predicted_sr = float(self.predicted_ensemble[0][0][1])
                predicted_mr = float(self.predicted_ensemble[0][0][0])

            elif self.ensemble_prediction_combination_mode in ["average", "voting"]:
                # Voting will also get a predicted_sr and mr calculated by averaging, since it can be useful to also
                #  use the averaged predictions in that case. Sometimes 4 models could yield outputs of a class as
                #  0.51, while a single model could classify the class as 0.1. The average will then be < 0.5, leading
                #  to a different predicted class based on averaging, compared to voting.
                predicted_sr = np.average([float(i[0][1]) for i in self.predicted_ensemble])
                predicted_mr = np.average([float(i[0][0]) for i in self.predicted_ensemble])

            if self.ensemble_prediction_combination_mode == "voting":
                # Index will be 1 for argmax in case of SR prediction, hence the sum of the argmaxes gives votes in
                #  favour of SR.
                # In case of having an even number of models and an equal number of votes for both classes, the class
                #  will be set as neither SR *nor* MR.

                votes = np.argmax([float(i[0]) for i in self.predicted_ensemble])
                pred_category_sr = 1 if np.sum(votes) > self.ensemble_model_count // 2 else 0
                pred_category_mr = 1 if np.sum(votes) < self.ensemble_model_count // 2 else 0

                reco_event_bs.Set("RingCountingVotingSRPrediction", pred_category_sr)
                reco_event_bs.Set("RingCountingVotingMRPrediction", pred_category_mr)
        else:
            predicted_sr = float(self.predicted[0][1])
            predicted_mr = float(self.predicted[0][0])

        reco_event_bs.Set("RingCountingSRPrediction", predicted_sr)
        reco_event_bs.Set("RingCountingMRPrediction", predicted_mr)

###################
# ↓ Boilerplate ↓ #
###################


thistool = RingCounting()


def SetToolChainVars(m_data_in, m_variables_in, m_log_in):
    return thistool.SetToolChainVars(m_data_in, m_variables_in, m_log_in)


def Initialise():
    return thistool.Initialise()


def Execute():
    return thistool.Execute()


def Finalise():
    return thistool.Finalise()
