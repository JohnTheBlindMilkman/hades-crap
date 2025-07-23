#ifndef FemtoMerger_hxx
    #define FemtoMerger_hxx

    #include <vector>
    #include <string>
    #include <iostream>
    #include "TFile.h"

    namespace Mixing
    {
        enum class Differenciate{No,Yes};

        struct Extent
        {
            std::string extension;
            std::vector<std::string> elements;
            Differenciate diff;
        };
        struct Histogram
        {
            std::string nameSignature;
            std::vector<Extent> extents;
        };
        namespace Detail
        {
            template <typename T>
            struct HistogramHandler
            {
                Histogram histogram;
                std::vector<std::pair<std::string,T*> > allHistograms;
            };
        }
        
        template <typename T>
        class FemtoMerger
        {
            private:
                std::string m_filePath;
                std::vector<Detail::HistogramHandler<T> > m_tasks;

                std::vector<std::string> CreateCombinations(const Histogram &hist)
                {
                    std::vector<std::string> allNames = {""};
                    for (const auto &extent : hist.extents)
                    {
                        std::vector<std::string> names;
                        for (const auto &name : allNames)
                        {
                            for (const auto &element : extent.elements)
                            {
                                names.push_back(name + element);
                            }
                        }
                        allNames = names;
                    }

                    for (const auto &name : allNames)
                        std::cout << name << "\t";

                    std::cout << "\n";
                    return allNames;
                }
                std::vector<std::pair<std::string,T*> > ReadData(const std::string &filePath, const Histogram &hist, const std::vector<std::string> &histExtents)
                {
                    std::vector<std::pair<std::string,T*> > allHistos;
                    TFile inpFile(filePath.c_str());

                    for (const auto &extent : histExtents)
                    {
                        allHistos.push_back(std::make_pair(extent,inpFile.Get<T>((hist.nameSignature + extent).c_str())));
                    }

                    for (auto hist : allHistos)
                        std::cout << hist.second->GetName() << "\t";

                    std::cout << "\n";
                    return allHistos;
                }

            public:
                FemtoMerger(const std::string &filePath, const std::vector<Histogram> &histograms) : m_filePath(filePath)
                {
                    for (const auto &hist : histograms)
                    {
                        m_tasks.push_back({hist,ReadData(filePath,hist,CreateCombinations(hist))});
                    }
                }

                std::vector<T*> MergeHistograms()
                {
                    std::vector<T*> mergedHistograms;
                    for (const auto &task : m_tasks)
                    {
                        for (std::size_t i = 0; i < task.histogram.extents.size(); ++i)
                        {
                            if (task.histogram.extents.at(i).diff == Differenciate::Yes)
                            {
                                std::vector<T*> histograms;
                                for (std::size_t j = 0; j < task.histogram.extents.at(i).elements.size(); ++j)
                                {
                                    std::vector<std::pair<std::string,T*> > histSubset;
                                    std::copy_if(m_tasks.at(i).allHistograms.begin(),
                                        m_tasks.at(i).allHistograms.end(),std::back_inserter(histSubset),
                                        [&](const std::pair<std::string,T*> &elem)
                                        {
                                            return (elem.first.at(i) == task.histogram.extents.at(i).elements.at(j));
                                        });

                                    T *tmpHist = nullptr;
                                    for (const auto &[name,hist] : histSubset)
                                    {
                                        if (hist != nullptr)
                                        {
                                            if (tmpHist == nullptr)
                                            {
                                                tmpHist = hist;
                                                tmpHist->SetName((task.histogram.nameSignature + task.histogram.extents.at(i).extension + std::to_string(j)).c_str());
                                            }
                                            else
                                            {
                                                tmpHist->Add(hist);
                                            }
                                        }
                                    }
                                    histograms.push_back(tmpHist);
                                }

                                std::move(histograms.begin(),histograms.end(),std::back_inserter(mergedHistograms));
                            }
                        }

                        T *histInteg = nullptr;
                        for (const auto &[name,hist] : task.allHistograms)
                        {
                            if (hist != nullptr)
                            {
                                if (histInteg == nullptr)
                                {
                                    histInteg = hist;
                                    histInteg->SetName((task.histogram.nameSignature + "Integ").c_str());
                                }
                                else
                                {
                                    histInteg->Add(hist);
                                }
                            }
                        }
                        mergedHistograms.push_back(histInteg);
                    }

                    return mergedHistograms;
                } // this is way too much nesting, but for now this will stay like that :/
        };
    } // namespace Mixing
    
#endif