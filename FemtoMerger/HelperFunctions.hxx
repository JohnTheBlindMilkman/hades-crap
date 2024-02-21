#ifndef HelperFunctions_hxx
    #define HelperFunctions_hxx

    #include <iostream>

    namespace JJFemtoMerger
    {
        class HelperFunctions
        {
            private:
                /* data */
            public:
                HelperFunctions() = default;
                ~HelperFunctions() = default;

                static void PrintHelp()
                {
                    std::cout << "Usage:\n";
                    std::cout << "./femtoMerge INPUTFILE SIGNATURE1 SIGNATURE2 KTMAX YMAX PSIMAX\n\n";
                    std::cout << "INPUTFILE - path to the analysis output .root file \n";
                    std::cout << "SIGNATURE1 - name used for signal histograms (without the numbers from FemtoMixer)\n";
                    std::cout << "SIGNATURE2 - name used for background histograms (without the numbers from FemtoMixer)\n";
                    std::cout << "KTMAX - number of the last k_T bin used in analysis and passed to FemtoMixer\n";
                    std::cout << "KTMAX - number of the last rapidity bin used in analysis and passed to FemtoMixer\n";
                    std::cout << "KTMAX - number of the last azimuthal angle w.r.t. the EP bin used in analysis and passed to FemtoMixer\n";
                };
        };
    }

#endif
