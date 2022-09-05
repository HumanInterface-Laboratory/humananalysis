from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import numpy as np
from scipy.signal import find_peaks
from scipy import interpolate
from scipy import integrate
import json
import os
from scipy.stats import linregress
from .py_ecg_detectors.ecgdetectors import Detectors
from .py_ecg_detectors.hrv import HRV
import scipy


class Biomarker():

    def __init__(self, DataPath: str, DeviceName: str) -> None:
        """データのパスからpandas形式で読み込み

        Args:
            DataPath (str): データファイルのパス
            DeviceName (str): 計測機器名[Nihonkoden, Nexus, Biolog]

        Examples:

        >>> Bio = Biomarker("./data.txt", DeviceName="Nexus")
        >>> Bio.showGraph()
        """

        if DeviceName == "Nihonkoden":
            with open(DataPath, 'r', encoding='cp932', newline='', errors="replace") as f:
                lines = list(f.read().splitlines())

                for i, line in enumerate(lines[:20]):
                    if line.startswith("Interval"):
                        if line.split('=')[-1] == "1.0ms":
                            self.Interval = 0.001  # 1ms
                        else:
                            print("Enter Interval time[ms]")
                            self.Interval = float(input())/1000

                    if line.startswith("#Address"):
                        self.Address = list(line.split('=')[-1].split(','))  # ex. [A1,A2,A3,A5,Event]
                        self.Address[-1] = "Event"
                    elif line.startswith("ChannelTitle"):
                        self.ChannelTitle = list(line.split('=')[-1].split(','))
                        self.ChannelTitle[-1] = "Event"
                        # print(self.ChannelTitle)

                    elif line.startswith("#Data="):
                        DataStartIdx = i
                        pass

            df = pd.read_csv(DataPath, header=None, sep=',', names=self.Address, encoding='cp932', skiprows=DataStartIdx+1, dtype="object")
            self.DataFrame = df

            self.EventsDF = df[df['Event'].str.startswith('#*', na=False)]['Event']

            self.MarkersDF = df[df['Event'].str.startswith('#* M', na=False)]['Event']
        elif DeviceName == "Nexus":
            with open(DataPath, 'r', encoding='cp932', newline='', errors="replace") as f:
                lines = list(f.read().splitlines())

                for i, line in enumerate(lines[:20]):
                    if line.startswith("Output rate:"):
                        if line.split('\t')[1] == "128":
                            self.Interval = 1/128  # 1ms
                        else:
                            print("Enter Output rate[Samples/sec]")
                            self.Interval = 1/int(input())

                    elif line.startswith("hh:mm:ss"):
                        self.Unit = list(line.split('\t'))  # ex. [hh:mm:ss	32 SPS	32 SPS	128 SPS	]

                    elif line.startswith("TIME"):
                        self.Address = list(line.split('=')[-1].split('\t'))
                        # del self.Address[-1]  # delete "Segments"
                        # print(self.ChannelTitle)

                    elif line.startswith("00:00:00"):
                        DataStartIdx = i-1
                        pass

            df = pd.read_table(DataPath, header=None, names=self.Address, encoding='cp932', skiprows=DataStartIdx+1, dtype="object")
            del df["Segments"]
            del df[""]
            del self.Address[-2:]  # delete "Segments"

            df.drop(df.tail(1).index, inplace=True)
            oldidx = 0
            oldtime = df["TIME"][0]
            timestamp = [0]
            for i, time in enumerate(df["TIME"]):
                if oldtime != time:
                    incremental = 1/(i - oldidx)
                    for j in range(oldidx, i-1):
                        timestamp.append(timestamp[-1] + incremental)
                    now = list(map(float, df["TIME"][i].split(":")))
                    timestamp.append(now[0]*3600 + now[1]*60 + now[2])
                    oldidx = i
                    oldtime = time
            incremental = 1/(i - oldidx)
            for j in range(oldidx, i):
                timestamp.append(timestamp[-1] + incremental)
            df["TIME"] = timestamp

            self.DataFrame = df
            self.EventsDF = df[~df['Events'].isna()]['Events']
            self.MarkersDF = df[~df['Events'].isna()]['Events']
        elif DeviceName == "Biolog":
            with open(DataPath, 'r', encoding='cp932', newline='', errors="replace") as f:
                lines = list(f.read().splitlines())
                lines[0].split(",")
                self.Address = lines[0].split(",")
                self.Unit = lines[2].split(",")
                DataStartIdx = 2

                self.Interval = (int(lines[4].split(",")[0].split(".")[-1]) - int(lines[3].split(",")[0].split(".")[-1]))//10

            df = pd.read_csv(DataPath, header=None, names=self.Address, encoding='cp932', skiprows=DataStartIdx+1, dtype="object")

            timestamp = []
            for i, time in enumerate(df["経過時間"]):
                now = list(map(int, time.split(".")))
                timestamp.append(now[0]*24*60*60*1000 + now[1]*60*60*1000 + now[2]*60*1000 + now[3]*1000 + now[4]/10)  # d.hh.mm.ss.msec
            df["経過時間"] = timestamp

            self.DataFrame = df
            # self.EventsDF = df[~df['Events'].isna()]['Events']
            # self.MarkersDF = df[~df['Events'].isna()]['Events']
            self.EventsDF = pd.DataFrame([])
            # TODO: DCH1以外に対応
            if "DCH1" in df.columns:
                self.MarkersDF = df[df['DCH1'].str.startswith('1', na=False)]['DCH1']
            else:
                self.MarkersDF = pd.DataFrame([])
        else:
            raise ValueError("Please enter the correct device name\n\nNihonkoden\nNexus\nBiolog")
        print("----------------------------------")
        print(f"{DeviceName=}")
        print("maker: index")
        print([[maker, index] for maker, index in zip(self.MarkersDF.values, self.MarkersDF.index)])
        # print(f"{self.EventsDF.values=}")
        # print(f"{self.EventsDF.index=}")
        print("----------------------------------")

    def calLFHF(
            self, column: str = "A5", starts: list = [],
            ends: list = [],
            prominence: float = 1, height: float = None, re_sampling_freq: float = 1,
            plot: bool = True, plot_vorbose: bool = False, revECG: bool = False) -> list:
        """自作のLF/HF算出関数(deprecated:非推奨)そのうち消す

        正しく算出できているか怪しいのでcalLFHF_v2()またはanalyzehrv(method="fanalysis")を使って


        """
        if len(starts) != len(ends):
            raise Exception('starts and ends are must be same length')
        LFHFlist = []
        for start, end in zip(starts, ends):
            ECGSignal = self.DataFrame[column].iloc[start:end].astype('float').values
            if revECG:
                ECGSignal = -ECGSignal
            # 以下岡野LFHF3を使用
            peak_index_array, _ = find_peaks(ECGSignal, prominence=prominence, height=height)
            time = np.arange(0, len(ECGSignal)*self.Interval, self.Interval)
            peak_time_array = peak_index_array*self.Interval

            if plot:
                plt.figure()
                plt.plot(time, ECGSignal, label="ECGSignal")
                plt.plot(peak_time_array, ECGSignal[peak_index_array], "ob", label="peaks")
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.title("peaks on ECGsignal")
                plt.tight_layout()
                # plt.xlim([1000,3000])
                plt.show()

            rr_intervals = (peak_time_array[1:]-peak_time_array[:-1])
            rr_times = peak_time_array[:-1]
            if plot_vorbose:
                plt.figure(figsize=(10, 5))
                plt.plot(rr_times, rr_intervals)
                plt.grid(True)
                # plt.ylim(0, )
                plt.show()

            # print(f"{rr_times=}")

            N = int(np.floor(rr_times[-1]-rr_times[0]/re_sampling_freq))  # re_sampling_freq で何点取れるか

            # print(N)

            f = interpolate.interp1d(rr_times, rr_intervals, kind='linear')

            t = np.arange(0, N, 1)
            y = f(rr_times[0]+t)
            if plot_vorbose:
                plt.figure()
                plt.plot(t, y)
                plt.grid(True)
                # plt.ylim(0, )
                plt.show()

            Y = np.fft.fft(y)
            fq = np.fft.fftfreq(len(y), d=1/re_sampling_freq)
            Y_abs_amp = abs(Y*Y)

            if plot_vorbose:
                plt.plot(fq[:N//2], Y_abs_amp[:N//2])
                plt.ylim(0, 8)
                plt.show()

            start_index = np.where(fq > 0.05)[0][0]
            end_index = np.where(fq > 0.15)[0][0]
            LF = integrate.cumtrapz(Y_abs_amp[start_index:end_index], fq[start_index:end_index])

            start_index = np.where(fq > 0.15)[0][0]
            end_index = np.where(fq > 0.40)[0][0]
            HF = integrate.cumtrapz(Y_abs_amp[start_index:end_index], fq[start_index:end_index])

            print("LF/HFは", LF[-1]/HF[-1])
            LFHFlist.append(LF[-1]/HF[-1])

        return LFHFlist

    def calLFHF_v2(
            self, column: str = "A5", starts: list = [],
            ends: list = [],
            prominence: float = 1, height: float = None, re_sampling_freq: float = 1,
            plot: bool = True, plot_vorbose: bool = False) -> list:

        if len(starts) != len(ends):
            raise Exception('starts and ends are must be same length')
        LFHFlist = []
        detectors = Detectors(int(1/self.Interval))
        hrv = HRV(int(1/self.Interval))
        for start, end in zip(starts, ends):
            ECGSignal = self.DataFrame[column].iloc[start:end].astype('float').values

            peak_index_array, filterd_signal = detectors.two_average_detector(ECGSignal)
            peak_index_array = np.array(peak_index_array)
            peak_index_array = peak_index_array[1:-1]
            time = np.arange(0, len(filterd_signal)*self.Interval, self.Interval)
            peak_time_array = peak_index_array*self.Interval

            if plot:
                plt.figure()
                plt.plot(time, filterd_signal, label="ECGSignal")
                plt.plot(peak_time_array, filterd_signal[peak_index_array], "ob", label="peaks")
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.title("peaks on ECGsignal")
                plt.tight_layout()
                # plt.xlim([1000,3000])
                plt.show()
            lfhf = hrv.fAnalysis(rr_samples=peak_index_array)

            print("LF/HFは", lfhf)
            LFHFlist.append(lfhf)

        return LFHFlist

    def analyzehrv(
            self, column: str = "A5", starts: list = [],
            ends: list = [],
            method: str = 'fAnalysis',
            detector_option: dict = {},
            plot: bool = True) -> list:
        """引数で指定された解析による結果を返す

        Args:
            column (str): 心電図のチャンネルにあたるカラム
            starts (List): 解析したい信号の開始インデックスのリスト
            ends (List): 解析したい信号の終了インデックスのリスト
            method (str): 心拍変動解析の手法名
            plot (bool): Rピークの検出結果のグラフを表示するかどうか

        Returns:
            List : 各start, end区間の指標のリスト

        Examples:

        >>> Bio = Biomarker("./data.txt", DeviceName="Nihonkoden")
        >>> Bio.analyzehrv("A5", [0], [-1], "fanalysis")
        """
        if len(starts) != len(ends):
            raise Exception('starts and ends are must be same length')

        hrv_index_list = []
        detectors = Detectors(int(1/self.Interval))
        hrv = HRV(int(1/self.Interval))

        for start, end in zip(starts, ends):
            ECGSignal = self.DataFrame[column].iloc[start:end].astype('float').values

            peak_index_array, filterd_signal = detectors.two_average_detector(ECGSignal, **detector_option)
            peak_index_array = np.array(peak_index_array)
            time = np.arange(0, len(filterd_signal)*self.Interval, self.Interval)
            peak_time_array = peak_index_array*self.Interval

            if plot:
                plt.figure()
                plt.plot(time, filterd_signal, label="ECGSignal")
                plt.plot(peak_time_array, filterd_signal[peak_index_array], "ob", label="peaks")
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.title("peaks on ECGsignal")
                plt.tight_layout()
                # plt.xlim([1000,3000])
                plt.show()

            # hrvの中のメソッドを取得、実行
            hrv_index = getattr(hrv, method)(rr_samples=peak_index_array)

            print(f"{method}は", hrv_index)
            hrv_index_list.append(hrv_index)

        return hrv_index_list

    def calHR(
            self, column: str = "A5", starts: list = [],
            ends: list = [],
            plot: bool = True, plot_vorbose: bool = False) -> list:
        """心拍数算出

        Args:
            column (str): 心電図のチャンネルにあたるカラム
            starts (List): 解析したい信号の開始インデックスのリスト
            ends (List): 解析したい信号の終了インデックスのリスト
            method (str): 心拍変動解析の手法名
            plot (bool): Rピークの検出結果のグラフを表示するかどうか

        Returns:
            List : 各start, end区間の心拍数のリスト

        Examples:

        >>> Bio = Biomarker("./data.txt", DeviceName="Nihonkoden")
        >>> print(Bio.calHR("A5", [0], [-1]))
        """
        if len(starts) != len(ends):
            raise Exception('starts and ends are must be same length')
        LFHFlist = []
        detectors = Detectors(int(1/self.Interval))
        hrv = HRV(int(1/self.Interval))
        for start, end in zip(starts, ends):
            ECGSignal = self.DataFrame[column].iloc[start:end].astype('float').values

            peak_index_array = detectors.two_average_detector(ECGSignal)
            peak_index_array = np.array(peak_index_array) - 50
            time = np.arange(0, len(ECGSignal)*self.Interval, self.Interval)
            peak_time_array = peak_index_array*self.Interval

            if plot:
                plt.figure()
                plt.plot(time, ECGSignal, label="ECGSignal")
                plt.plot(peak_time_array, ECGSignal[peak_index_array], "ob", label="peaks")
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.title("peaks on ECGsignal")
                plt.tight_layout()
                # plt.xlim([1000,3000])
                plt.show()
            HR = np.mean(hrv.HR(rr_samples=peak_index_array))

            print("HRは", HR)
            LFHFlist.append(HR)

        return LFHFlist

    def showGraph(self, columns: list = [], setplot: bool = True, divplot: bool = True):
        """データの可視化

        最後の列はマーカーが入っている想定なので無視していることに注意

        Args:
            columns (List): グラフを表示するチャンネル名
            setplot (bool): 一つのウィンドウにまとめて表示
            divplot (bool): 個別にウィンドウに表示


        Examples:

        >>> Bio = Biomarker("./data.txt", DeviceName="Nihonkoden")
        >>> Bio.showGraph()
        """
        df = self.DataFrame
        if columns == []:
            columns = self.Address[:-1]
        if setplot:
            _, axes = plt.subplots(len(columns), 1)
            for column, ax in zip(columns, axes):
                ax.plot(df.index * self.Interval, df[column].astype('float'), zorder=1)  # 筋電データ
                for idx in self.MarkersDF.index:
                    ax.axvline(idx * self.Interval, ls="-", color="red")
                ax.set_xlabel("time [sec]")
                ax.set_ylabel(f"{column}")

        if divplot:
            for column in columns:
                plt.figure()
                plt.plot(df.index * self.Interval, df[column].astype('float'), zorder=1)  # 筋電データ
                for idx in self.MarkersDF.index:
                    plt.axvline(idx * self.Interval, ls="-", color="red")
                plt.xlabel("time [sec]")
                plt.ylabel(f"{column}")

        plt.show()

        pass

    def getDataFrame(self, columns: list = []):
        df = self.DataFrame
        if columns == []:
            columns = self.Address[:-1]
        return df.loc[:, columns]

    def getRowSignal(self, column="", starts: list = [], ends: list = []) -> list:
        df = self.DataFrame
        signals = []

        for start, end in zip(starts, ends):
            signal = df[column].iloc[start:end].astype(float).values.flatten()
            signals.append(signal)
        return signals

    def calRootMeanSquareSignal(self, columns=[], starts: list = [], ends: list = [], window_time_ms: int = 50, plot: bool = False) -> list:
        def _rms(d): return np.sqrt((d ** 2).sum() / d.size)
        df = self.DataFrame
        signals = []
        for column in columns:
            for start, end in zip(starts, ends):
                window_size = int(window_time_ms/(self.Interval*1000))
                # df["rms_signal"] = df[column].iloc[start:end].rolling(window=window_size, min_periods=1, center=False).apply(_rms)
                # rms_signal = df["rms_signal"].values
                rms_signal = df[column].iloc[start:end].rolling(window=window_size, min_periods=1, center=False).apply(_rms).values
                signals.append(rms_signal)
                if plot:
                    # print(len(rms_signal), start, end)
                    self.__dfplot(rms_signal, start, end)
        return signals

    def calMVC(self, columns=[], starts: list = [], ends: list = [], window_time_ms: int = 50, plot: bool = False) -> list:
        def _rms(d): return np.sqrt((d ** 2).sum() / d.size)
        df = self.DataFrame
        signals = []
        for column in columns:
            for start, end in zip(starts, ends):
                window_size = int(window_time_ms/(self.Interval*1000))
                df["rms_signal"] = df[column].iloc[start:end].rolling(window=window_size, min_periods=1, center=False).apply(_rms)
                rms_signal = df["rms_signal"].dropna().values
                # signals.append(rms_signal)
                signal = df["rms_signal"].iloc[start:end].rolling(window=3000, min_periods=1, center=False).apply(sum).values
                peak_index_array, _ = find_peaks(signal, prominence=30, height=None, distance=None)
                peak_index_array_r, _ = find_peaks(-signal, prominence=30, height=None, distance=None)
                plt.plot(range(len(signal)), signal)
                plt.plot(peak_index_array, signal[peak_index_array], "ob")
                plt.plot(peak_index_array_r, signal[peak_index_array_r], "ob")
                plt.show()

                # rms_time = np.arange(0, len(rms_signal)*self.Interval, self.Interval)
                # time = np.arange(0, len(signal)*self.Interval, self.Interval)
                # peak_time_array = peak_index_array*self.Interval
                # peak_time_array_r = peak_index_array_r*self.Interval
                starts = peak_index_array[0]
                print(f"{peak_index_array=}")
                ends = peak_index_array[0]+3000
                print(starts, ends)
                signals.append(rms_signal[starts-3000:ends-3000])
                if plot:
                    plt.figure()
                    # print(list(range(len(rms_signal)))[0:50])
                    print(rms_signal)
                    plt.plot(list(range(len(rms_signal))), rms_signal, label="-signal")
                    plt.plot(peak_index_array-3000, rms_signal[peak_index_array-3000], "ob", label="peaks")
                    plt.plot(peak_index_array_r-3000, rms_signal[peak_index_array_r-3000], "ob", label="peaks")
                    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                    plt.title("peaks on signal")
                    plt.tight_layout()
                    plt.show()
                    plt.clf()
                    plt.close()

        return signals

    # def calEMGActiveStartEndTime(self, signal=[], prominence: float = 1, height: float = None, distance: float = 20, plot: bool = False) -> list:
    #     starts = []
    #     ends = []

    #     signal = np.array(signal)
    #     peak_index_array, _ = find_peaks(-signal, prominence=prominence, height=height, distance=distance)
    #     time = np.arange(0, len(signal)*self.Interval, self.Interval)
    #     peak_time_array = peak_index_array*self.Interval
    #     starts.append([start for start in peak_index_array[:-1]])
    #     ends.append([end for end in peak_index_array[1:]])

    #     if plot:
    #         plt.plot(time, signal, label="-signal")
    #         plt.plot(peak_time_array, signal[peak_index_array], "ob", label="peaks")
    #         plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    #         plt.title("peaks on signal")
    #         plt.tight_layout()
    #         plt.show()
    #     return starts, ends

    def calAverageAmplitude(self, signal=[], starts=[0], ends=[-1], plot: bool = False):

        average_amplitudes = []
        for start, end in zip(starts, ends):
            # print(start, end)
            _signal = signal[start:end]
            # print(signal)
            # print(signal)
            if plot:
                plt.plot(_signal)
                plt.show()
            average_amplitude = sum(_signal)/len(_signal)
            average_amplitudes.append(average_amplitude)

        return average_amplitudes

    def calNumOfPeaks(
            self, signals=[],
            MVC: float = None, prominence: float = 1, height: float = None, distance: float = 20, plot: bool = False, threshold: float = None) -> list:
        num_of_peaks = []
        for signal in signals:
            peak_index_array, _ = find_peaks(signal, prominence=prominence, height=height, distance=distance)
            if plot:
                plt.plot(signal, label="-signal")
                plt.plot(peak_index_array, signal[peak_index_array], "ob", label="peaks")
                if MVC is not None:
                    plt.axhline(MVC)
                plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
                plt.title("peaks on signal")
                plt.tight_layout()
                plt.show()
            if threshold is not None:
                peak_index_array = [peak for peak in peak_index_array if peak >= threshold]
            num_of_peaks.append(len(peak_index_array))
        return num_of_peaks

    def annotation(self, columns: list = [], starts: list = [], ends: list = [], filename: str = ""):
        def select_callback(eclick, erelease):
            """
            Callback for line selection.
            *eclick* and *erelease* are the press and release events.
            """
            x1, y1 = eclick.xdata, eclick.ydata
            x2, y2 = erelease.xdata, erelease.ydata
            # print(f"({x1:3.2f}, {y1:3.2f}) --> ({x2:3.2f}, {y2:3.2f})")
            # print(f"The buttons you used were: {eclick.button} {erelease.button}")

        def toggle_selector(event):
            print('Key pressed.')
            # if event.key == 't':
            #     name = type(selector).__name__
            #     if selector.active:
            #         print(f'{name} deactivated.')
            #         selector.set_active(False)
            #     else:
            #         print(f'{name} activated.')
            #         selector.set_active(True)
            if event.key == "a":
                xmin, xmax, ymin, ymax = selector.extents
                xmin, xmax, ymin, ymax = int(xmin), int(xmax), int(ymin), int(ymax)
                ax.axvline(xmin, color="r")
                ax.axvline(xmax, color="pink")
                print(xmin, xmax)
                if xmin not in dic[column]["chewing_starts"]:
                    dic[column]["chewing_starts"].append(xmin)
                    dic[column]["chewing_ends"].append(xmax)
                with open(f"./addData/{filename}.json", "w") as f:
                    json.dump(dic, f, indent=1)
                plt.draw()
            if event.key == "t":
                xmin, xmax, ymin, ymax = selector.extents
                xmin, xmax, ymin, ymax = int(xmin), int(xmax), int(ymin), int(ymax)
                ax.axvline(xmin, color="gold")
                ax.axvline(xmax, color="lemonchiffon")
                print(xmin, xmax)
                if xmin not in dic[column]["MVC_starts"]:
                    dic[column]["MVC_starts"].append(xmin)
                    dic[column]["MVC_ends"].append(xmax)
                with open(f"./addData/{filename}.json", "w") as f:
                    json.dump(dic, f, indent=1)
                plt.draw()
        df = self.DataFrame

        if os.path.exists(f"./addData/{filename}.json"):
            with open(f"./addData/{filename}.json", "r") as f:
                try:
                    dic = json.load(f)
                    for column in columns:
                        try:
                            dic[column]
                        except:
                            dic[column] = {"chewing_starts": [], "chewing_ends": [], "MVC_starts": [], "MVC_ends": []}
                except:
                    dic = {column: {"chewing_starts": [], "chewing_ends": [], "MVC_starts": [], "MVC_ends": []} for column in columns}
        else:
            dic = {column: {"chewing_starts": [], "chewing_ends": [], "MVC_starts": [], "MVC_ends": []} for column in columns}
            Path(f"./addData/{filename}.json").touch(exist_ok=False)

        for column in columns:
            for start, end in zip(starts, ends):
                signal = df[column].iloc[start:end].astype(float).values
                fig = plt.figure(constrained_layout=True)
                ax = fig.subplots()
                ax.plot(range(len(signal)), signal)
                ax.set_title(filename+column+"\nPress 't' to mark MVC"+"\nPress 'a' to mark sosyaku")

                for idx in df.iloc[start:end][df.iloc[start:end]['Event'].str.startswith('#* M', na=False)].index:
                    plt.axvline(idx, ls="-", color="green")
                # print(dic)
                chewing_starts, chewing_ends = dic[column]["chewing_starts"], dic[column]["chewing_ends"]
                MVC_starts, MVC_ends = dic[column]["MVC_starts"], dic[column]["MVC_ends"]
                for _start in chewing_starts:
                    plt.axvline(_start, ls="-", color="gray")
                for _end in chewing_ends:
                    plt.axvline(_end, ls="-", color="lightgray")
                for _start in MVC_starts:
                    plt.axvline(_start, ls="-", color="gray")
                for _end in MVC_ends:
                    plt.axvline(_end, ls="-", color="lightgray")

                selector = RectangleSelector(ax, select_callback, useblit=True, button=[1, 3],  # disable middle button
                                             minspanx=5, minspany=5,
                                             spancoords='pixels',
                                             interactive=True)
                fig.canvas.mpl_connect('key_press_event', toggle_selector)

                plt.show()

    def detect_emg_activations(self, emg_signal, sample_rate, smooth_level=20, threshold_level=10,
                               time_units=False, volts=False, resolution=None, device="biosignalsplux",
                               plot_result=False):
        """
        -----
        Brief
        -----
        Python implementation of Burst detection algorithm using Teager Kaiser Energy Operator.
        -----------
        Description
        -----------
        Activation events in EMG readings correspond to an increase of muscular activity, namely, from inaction to action.
        These events are characterised by an increase in electric potential that returns to the initial values when the
        muscle returns to a state of inaction.
        This function detects activation events using the Teager Kaiser Energy Operator.
        ----------
        Parameters
        ----------
        emg_signal : list
            List of EMG acquired samples.
        sample_rate : int
            Sampling frequency.
        smooth_level : number
            Defines a percentage proportional to the smoothing level, i.e. the bigger this value is,
            the more smoothed is the signal.
        threshold_level : number
            Specification of the single threshold position, used for distinguishing between activation
            (above) and inactivation samples (below).
        time_units : boolean
            If True this function will return the Burst begin and end positions in seconds.
        volts : boolean
            If True, then the conversion of raw units to mV will be done. Resolution need to be
            specified.
        resolution : int
            Selected resolution for data acquisition.
        device : str
            Specification of the device category.
        plot_result : boolean
            If True it will be presented a graphical representation of the detected burst in the EMG
            signal.
        Returns
        -------
        out : bursts begin (ndarray), bursts end (ndarray)
            Begin and end of bursts (sample number or time instant in seconds).
        smooth_signal: list
            It is returned the smoothed EMG signal (after the processing steps intended to simplify the
            signal).
        threshold_level: float
            The value of the detection threshold used to locate the begin and end of each muscular
            activation period.
        """

        if volts is True:
            if resolution is not None:
                emg_signal = raw_to_phy("EMG", device, emg_signal, resolution, option="mV")
                units = "mV"
            else:
                raise RuntimeError(
                    "For converting raw units to mV is mandatory the specification of acquisition "
                    "resolution.")
        else:
            units = "Input Units"

        if time_units is True:
            time_units_str = "Time (s)"
            time = np.linspace(0, len(emg_signal) / sample_rate, len(emg_signal))
        else:
            time = np.linspace(0, len(emg_signal) - 1, len(emg_signal))
            time_units_str = "Sample Number"

        # ----------------------------------- Baseline Removal -----------------------------------------
        pre_pro_signal = np.array(emg_signal) - np.average(emg_signal)

        # ------------------------------------ Signal Filtering ----------------------------------------
        # low_cutoff = 10  # Hz
        # high_cutoff = 300  # Hz

        # # Application of the signal to the filter.
        # pre_pro_signal = _butter_bandpass_filter(pre_pro_signal, low_cutoff, high_cutoff, sample_rate)

        # ------------------------------ Application of TKEO Operator ----------------------------------
        tkeo = []
        for i, signal_sample in enumerate(pre_pro_signal):
            if i in (0, len(pre_pro_signal) - 1):
                tkeo.append(signal_sample)
            else:
                tkeo.append(np.power(signal_sample, 2) - (pre_pro_signal[i + 1] *
                                                          pre_pro_signal[i - 1]))

        # Smoothing level - Size of sliding window used during the moving average process (a function
        # of sampling frequency)
        smoothing_level = int((smooth_level / 100) * sample_rate)

        # --------------------------------- Signal Rectification ---------------------------------------
        rect_signal = np.absolute(tkeo)

        # ------------------------------ First Moving Average Filter -----------------------------------
        rect_signal = self.__moving_average(rect_signal, sample_rate / 10)

        # -------------------------------- Second Smoothing Phase --------------------------------------
        smooth_signal = []
        for i in range(0, len(rect_signal)):
            if smoothing_level < i < len(rect_signal) - smoothing_level:
                smooth_signal.append(np.mean(rect_signal[i - smoothing_level:i + smoothing_level]))
            else:
                smooth_signal.append(0)

        # ----------------------------------- Threshold -----------------------------------------------
        avg_pre_pro_signal = np.average(pre_pro_signal)
        std_pre_pro_signal = np.std(pre_pro_signal)

        threshold_level = avg_pre_pro_signal + self.__thres_norm_reg(threshold_level, smooth_signal,
                                                                     pre_pro_signal) * std_pre_pro_signal

        # Generation of a square wave reflecting the activation and inactivation periods.
        binary_signal = []
        for i in range(0, len(time)):
            if smooth_signal[i] >= threshold_level:
                binary_signal.append(1)
            else:
                binary_signal.append(0)

        # ------------------------------ Begin and End of Bursts --------------------------------------
        diff_signal = np.diff(binary_signal)
        act_begin = np.where(diff_signal == 1)[0]
        act_end = np.where(diff_signal == -1)[0]

        if time_units is True:
            time_begin = np.array(time)[act_begin]
            time_end = np.array(time)[act_end]
        else:
            time_begin = act_begin
            time_end = act_end

        # If plot is invoked by plot_result flag, then a graphical representation of the R peaks is
        # presented to the user.
        # if plot_result is True:
        #     plot([list(time), list(time)], [list(emg_signal), list(np.array(binary_signal) *
        #                                                            np.max(emg_signal))],
        #          yAxisLabel=["Data Samples (" + units + ")"] * 2,
        #          x_axis_label=time_units_str, legend_label=["EMG Signal", "Activation Signal"])

        return time_begin, time_end, smooth_signal, threshold_level

    def emg_fanalysis(self, signal=[], plot: bool = False):
        f, Pxx = scipy.signal.welch(signal, fs=1/self.Interval)
        # print(f, t, Zxx)
        # print(len(Pxx))
        if plot:

            plt.bar(f, Pxx)
            # plt.title(filename + " " + str(i))
            # plt.pcolormesh(t, f, np.abs(Zxx))  # , shading='gouraud'
            # plt.plot(f)
            plt.tight_layout()
            plt.show()

        LFB = 0
        MFB = 0
        HFB = 0
        # 15~45Hz,45~80Hz,80Hz~500Hz
        for _f, _Pxx in zip(f, Pxx):
            if 15 <= _f <= 45:
                LFB += _Pxx
            elif 45 < _f <= 80:
                MFB += _Pxx
            elif 80 < _f <= 500:
                HFB += _Pxx
            else:
                pass
        half_of_sumPxx = sum(Pxx)/2
        p = 0
        for i, _p in enumerate(Pxx):
            p += _p
            if half_of_sumPxx <= p:
                break
        IntermediateFrequency = f[i]
        _sum = sum((LFB, MFB, HFB))
        return LFB, MFB, HFB, LFB/_sum, MFB/_sum, HFB/_sum, IntermediateFrequency

    def __dfplot(self, signal, start, end):
        df = self.DataFrame
        plt.figure()
        plt.plot(df.iloc[start:end].index * self.Interval, signal, zorder=1)
        for idx in df.iloc[start:end][df.iloc[start:end]['Event'].str.startswith('#* M', na=False)].index:
            plt.axvline(idx * self.Interval, ls="-", color="red")
        plt.xlabel("time [msec]")
        plt.ylabel("amplified EMG [mV]")
        plt.show()

    def __moving_average(self, data, wind_size=3):
        """
        -----
        Brief
        -----
        Application of a moving average filter for signal smoothing.
        -----------
        Description
        -----------
        In certain situations it will be interesting to simplify a signal, particularly in cases where
        some events with a random nature take place (the random nature of EMG activation periods is
        a good example).
        One possible simplification procedure consists in smoothing the signal in order to obtain
        only an "envelope". With this methodology the analysis is mainly centered on seeing patterns
        in data and excluding noise or rapid events [1].
        The simplification can be achieved by segmenting the time series in multiple windows and
        from each window an average value of all the samples that it contains will be determined
        (dividing the sum of all sample values by the window size).
        A quick and efficient implementation (chosen in biosignalsnotebooks package) of the moving window
        methodology is through a cumulative sum array.
        [1] https://en.wikipedia.org/wiki/Smoothing
        ---------
        Parameters
        ----------
        data : list
            List of signal samples.
        wind_size : int
            Number of samples inside the moving average window (a bigger value implies a smoother
            output signal).
        Returns
        -------
        out : numpy array
            Array that contains the samples of the smoothed signal.
        """

        wind_size = int(wind_size)
        ret = np.cumsum(data, dtype=float)
        ret[wind_size:] = ret[wind_size:] - ret[:-wind_size]
        return np.concatenate((np.zeros(wind_size - 1), ret[wind_size - 1:] / wind_size))

    def __thres_norm_reg(self, threshold_level, signal, pre_smooth_signal):
        """
        Regression function that with a percent input gives an absolute value of the threshold
        level (used in the muscular activation detection algorithm).
        Converts a relative threshold level to an absolute value.
        ----------
        Parameters
        ----------
        threshold_level : int
            Percentage value that defines the absolute threshold level relatively to the maximum value
            of signal.
        signal : list
            List of EMG smoothed signal samples.
        pre_smooth_signal : list
            Original EMG samples.
        Returns
        -------
        out : float
            Threshold level in absolute format.
        """
        avg_signal = np.average(pre_smooth_signal)
        std_signal = np.std(pre_smooth_signal)

        threshold_0_perc_level = (-avg_signal) / float(std_signal)
        threshold_100_perc_level = (np.max(signal) - avg_signal) / float(std_signal)

        slope, b_coeff = linregress([0, 100], [threshold_0_perc_level, threshold_100_perc_level])[:2]
        return slope * threshold_level + b_coeff


def main():

    pass


if __name__ == "__main__":

    main()
