"""
Creates the test networks needed to reproduce our papers' results
"""
function create_enwl_4w_data()
    create_loads_files()
    create_lines_files()
    create_master_files()
end
"""
Creates updated load files for both test networks
"""
function create_loads_files()
# Open the input text file
    input_files  = [joinpath(@__DIR__, "network_data", "30load-feeder", "Loads.txt"),  joinpath(@__DIR__, "network_data", "eulvtf", "Loads.txt")]
    output_files = [joinpath(@__DIR__, "network_data", "30load-feeder", "Loads2.txt"), joinpath(@__DIR__, "network_data", "eulvtf", "Loads2.txt")]

    for (input_file, output_file) in zip(input_files, output_files)
        # Read the input file line by line
        lines = readlines(input_file)

        # Open the output file for writing
        open(output_file, "w") do io
            # Write header
            write(io, "! Imported from Loads.txt, with modified neutral terminal\n")
            write(io, "! ---------------------------------------------------------------------------\n")
            
            # Process each line and write the modified line to the output file
            for line in lines
                modified_line = modify_line_loads(line)
                write(io, modified_line)
            end
            write(io, "! ---------------------------------------------------------------------------\n")
        end
    end
end
"""
Changes the line entries of each load text file to accommodate the neutral conductor
"""
function modify_line_loads(line::String)
    parts = split(line)
    bus_part = parts[4] * ".4"
    # kW_value = parse(Float64, split(parts[6], "=")[2])
    # kW_new = kW_value * rand() * 3 # Adjust the multiplication factor as needed    
    return "New Load.$(parts[2][6:end]) Phases=1 Bus1=$bus_part kV=0.23 kW=1  PF=0.95\n"
end
"""
Changes the line entries of the lines.txt
"""
function modify_line_lines(line::String, is_first_line::Bool, output_filename::String)
    parts = split(line)
    
    # Extract line number from the format "New Line.LINE1"
    line_num_str = match(r"LINE(\d+)", parts[2]).captures[1]
    line_num = parse(Int, line_num_str)
        
    if occursin("eulvtf", output_filename)  
        linecode = occursin("_ug", output_filename) ? "uglv_120cu_xlpe/nyl/pvc_ug_4w_bundled" : "pluto"
    else # 30load network
        linecode = occursin("_ug", output_filename) ? "uglv_240al_xlpe/nyl/pvc_ug_4w_bundled" : "pluto"
    end
    
    if line_num == 328 && occursin("30load", output_filename) && occursin("_ug", output_filename)
        bus1 = "322.1.2.3.4"
    else
        bus1 = "$(split(parts[3], "=")[end]).1.2.3.4"    
    end
    bus2 = "$(split(parts[4], "=")[end]).1.2.3.4"    
    
    length_value = parse(Float64, split(parts[7], "=")[2])
    
    if is_first_line
        return """New Line.LINE0 Bus1=sourcebus.1.2.3.4 Bus2=1.1.2.3.4 phases=4 Linecode=$linecode Length=1 Units=m
                  New Line.$(parts[2][6:end]) Bus1=$bus1 Bus2=$bus2 phases=4 Linecode=$linecode Length=$length_value Units=m\n"""
    else
        return "New Line.$(parts[2][6:end]) Bus1=$bus1 Bus2=$bus2 phases=4 Linecode=$linecode Length=$length_value Units=m\n"
    end
end
"""
Creates `Lines_oh.txt` and `Lines_ug.txt` for both test cases
"""
function create_lines_files()
    input_files  = [joinpath(@__DIR__, "network_data", "30load-feeder", "Lines.txt"),  joinpath(@__DIR__, "network_data", "eulvtf", "Lines.txt")]

    for case in ["_ug", "_oh"]
        output_files = [joinpath(@__DIR__, "network_data", "30load-feeder", "Lines$case.txt"), joinpath(@__DIR__, "network_data", "eulvtf", "Lines$case.txt")]
        for (input_file, output_file) in zip(input_files, output_files)
            # Read the input file line by line
            lines = readlines(input_file)

            # Open the output file for writing
            open(output_file, "w") do io
                write(io, "! Imported from Lines.txt, with modified linecodes\n")
                if occursin("30load", output_file) && occursin("_ug", output_file) 
                    write(io, "! WARNING! change w.r.t. original ENWL data: line328 is from bus 322 to 329 instead of 326 to 329!!!\n")
                end
                write(io, "! ---------------------------------------------------------------------------\n")
                # Process each line and write the modified line to the output file
                for (index, line) in enumerate(lines)
                    modified_line = modify_line_lines(line, index == 1, output_file)
                    write(io, modified_line)
                end
                write(io, "! ---------------------------------------------------------------------------\n")
            end
        end
    end
end
"""
Creates master files for all four cases
"""
function create_master_files()
    output_files = [joinpath(@__DIR__, "network_data", "30load-feeder", "Master_oh.dss"), joinpath(@__DIR__, "network_data", "30load-feeder", "Master_ug.dss"),
                    joinpath(@__DIR__, "network_data", "eulvtf", "Master_oh.dss"), joinpath(@__DIR__, "network_data", "eulvtf", "Master_ug.dss")]
    
    common_header = """! Original network and loads data from ENWL - Low Voltage Network Solutions project (in OpenDSS format)
                      !   https://www.enwl.co.uk/go-net-zero/innovation/smaller-projects/low-carbon-networks-fund/low-voltage-network-solutions/ 
                      Clear
                      Set DefaultBaseFreq=50\n"""
    
    common_center = """New Reactor.Grounding phases=1 bus1=sourcebus.4 bus2=sourcebus.0 R=0.0 X=1E-10
                     Redirect ../../linecode_library/LineCode_impedances_reordered.dss
                     Redirect Loads.txt\n"""
    
    common_ending = """New Energymeter.substation Element=Line.LINE1 1
                     Set mode=Snap
                     Solve
                     Closedi"""

    for output_file in output_files
        network = occursin("eulvtf", output_file) ? "New Circuit.ENWL_network_1_Feeder_1_4wire BasekV=0.416 pu=1.00 ISC3=100000 ISC1=100000\n" : "New Circuit.ENWL_network_5_Feeder_4_4wire BasekV=0.416 pu=1.00 ISC3=100000 ISC1=100000\n"
        ug_or_oh = occursin("_oh", output_file) ? "Redirect Lines_oh.txt\n" : "Redirect Lines_ug.txt\n"
        open(output_file, "w") do io
            write(io, common_header)
            write(io, network)
            write(io, common_center)
            write(io, ug_or_oh)
            write(io, common_ending)
        end
    end
end

