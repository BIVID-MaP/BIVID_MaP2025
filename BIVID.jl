using Pkg
using CSV
using DataFrames
using BioAlignments
using BioSequences
using Statistics
using StringDistances
using Base.Threads
using Distributed
using Bio
using Bio.Seq
using ArgParse
using JSON
using BioAlignments
using BioSequences
using BioTools
using SAMTools
using Dates
function parse_args()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--fasta_path"
            help = "FASTA file path including both reference sequence and mutant sequences"
            default = "./Demo/Input_file/test_FASTA_G4I8.txt"
        "--sam_dir"
            help = "SAM file dir path mapped to reference"
            default = "./Demo/Input_file/input_sam"
    end
    return ArgParse.parse_args(s)  
end
 
function extract_reads_by_fasta_id(input_samfile::AbstractString, output_samfile::AbstractString, fasta_id::AbstractString)
    open(input_samfile, "r") do input_file
        open(output_samfile, "w") do output_file
            for line in eachline(input_file)
                # ヘッダー行はそのまま新しいファイルに書き込む
                if startswith(line, '@')
                    write(output_file, line, "\n")
                    continue
                end
                
                fields = split(line, '\t')
                read_id = fields[3]

                # リードのFasta IDを確認       
                if read_id == fasta_id
                    
                    write(output_file, line, "\n")
                end
            end
        end
    end
end
function filter_sam_by_base(samfile::String, n::Int, base::Char, outputfile::String)###n番目の文字にbaseの塩基が入る場合

    open(samfile, "r") do input
        open(outputfile, "w") do output
            # ヘッダーの書き込み
            while true 
                line = readline(input)
           
                if line[1] == '@'
                    write(output, line)
                    write(output, '\n')
                else
                    break
                end
            end
            # リードのフィルタリングと書き込み
            while !eof(input)###リード0のファイルではない  
                line = readline(input)
                fields = split(line, '\t')
   
                query = fields[10]  

                MAPPOS=parse(Int,fields[4])-1
 
                cigar= fields[6]
               
                if occursin("H",cigar)==true

                else
                    c_s=makeCIGARstring(cigar)
                    if length(c_s) > n-1-MAPPOS && MAPPOS < n 
                        
                        if n==1 && c_s[n-MAPPOS] == 'M' && query[n-MAPPOS] == base
               
                                write(output, line)
                                write(output, '\n')
            
                        elseif count_M(c_s[1:n-MAPPOS])==length(c_s[1:n-MAPPOS]) && query[n-MAPPOS] == base && c_s[n-MAPPOS]=='M'
              
                                write(output, line)
                                write(output, '\n')
                            
                        elseif count_D(c_s[1:n-1-MAPPOS]) >0 || count_I(c_s[1:n-1-MAPPOS]) >0
                        
           
                            search_n=n-count_D(c_s[1:n-1-MAPPOS])+count_I(c_s[1:n-1-MAPPOS])-MAPPOS
         
                                if query[search_n] == base && c_s[search_n]=='M'
                                    write(output, line)
                                    write(output, '\n')
                                end
                        end
                    end

                end

            end
        end
    end
end
# 使用例



function extract_reads_with_deletions(input_samfile::AbstractString, output_samfile::AbstractString, deletion_position_list)
    open(input_samfile, "r") do input_file
        open(output_samfile, "w") do output_file
            for line in eachline(input_file)
                if startswith(line, '@')
                    write(output_file, line, "\n")
                    continue
                end
                fields = split(line, '\t')
                cigar = fields[6]
                MAPPOS=parse(Int,fields[4])-1
                deletion_position_list2=map(x->x-MAPPOS,deletion_position_list)
                seq = fields[10]
                cigar_string=makeCIGARstring(cigar)
                if occursin("D",cigar_string)==true
                    count=0
                    for c in cigar_string
                        count+=1
                    
                        for delpos in deletion_position_list2
                            if c=='D' && count==delpos
                                write(output_file, line, "\n")
                            end
                        end
                    end
                end
            end
        end
    end
end

function extract_reads_with_alldeletions(input_samfile::AbstractString, output_samfile::AbstractString)
open(input_samfile, "r") do input_file
    open(output_samfile, "w") do output_file
        for line in eachline(input_file)
            if startswith(line, '@')
                write(output_file, line, "\n")
                continue
            end
            fields = split(line, '\t')
            cigar = fields[6]
            MAPPOS=parse(Int,fields[4])-1

            seq = fields[10]
            cigar_string=makeCIGARstring(cigar)

            # 欠損があるかどうかを確認
            if occursin("D",cigar_string)==true
                write(output_file, line, "\n")
            end
        end
    end
end
end

find_snv(ref::AbstractString, mut::AbstractString) = let
    @assert length(ref) == length(mut)
    i = findfirst(i -> ref[i] ≠ mut[i], 1:length(ref))
    i === nothing ? nothing : (i, ref[i], mut[i])
end


function read_fasta(path::AbstractString)
    seqs = Dict{String,String}()
    current = ""
    for line in eachline(path)
        if startswith(line, '>')
            current = strip(line[2:end])
            seqs[current] = ""
        else
            seqs[current] *= strip(line)
        end
    end
    return seqs
end


function make_divided_sam(fasta_path,sam_path,output_folder)

    seqs = read_fasta(fasta_path)

    ref_name = only(filter(name -> endswith(name, "_Ref"), keys(seqs)))
    ref_seq  = seqs[ref_name]
    mut_names = filter(name -> endswith(name, "_Ref") == false, keys(seqs))

    mut_pos_dict = Dict{String, Tuple{Int,Char,Char}}()
    for mut in mut_names
        snv = find_snv(ref_seq, seqs[mut])
        snv === nothing && continue
        idx, r, m = snv

        mut_pos_dict[mut] = (idx, r, m)
    end
    parent_samfile_identifier=split(split(sam_path,"/")[end],".sam")[1]
    seq_name=split(ref_name,"_")[1]
    target_sam="$seq_name.$parent_samfile_identifier.sam"
    extract_reads_by_fasta_id(sam_path, "$output_folder/$target_sam", ref_name)


    for (mutname, mut_values) in zip(keys(mut_pos_dict), values(mut_pos_dict))
        output_divided_sam=replace(target_sam,".sam"=>".$mutname.divided.sam")
        filter_sam_by_base("$output_folder/$target_sam", 
        mut_values[1],
        mut_values[3]
        , "$output_folder/$output_divided_sam")
    end

    current_sam="$output_folder/$target_sam"
    for (mutname, mut_values) in zip(keys(mut_pos_dict), values(mut_pos_dict))
        middle_reference_divided_sam=replace(target_sam,".sam"=>".$ref_name.$mutname.divided.middle.sam")
        filter_sam_by_base(current_sam, 
        mut_values[1],
        mut_values[2]
        , "$output_folder/$middle_reference_divided_sam")
        if current_sam != "$output_folder/$target_sam"
            rm(current_sam)
        end
        current_sam="$output_folder/$middle_reference_divided_sam"
    end
    output_reference_divided_sam=replace(target_sam,".sam"=>".$ref_name.divided.sam")
    final_file = "$output_folder/$output_reference_divided_sam"
    mv(current_sam, final_file,force=true)

    all_deletion_outputsam=replace(target_sam,".sam"=>".all_deletion.sam")
    variantpos_deletion_outputsam=replace(target_sam,".sam"=>".variantpos_deletion.sam")
    extract_reads_with_alldeletions("$output_folder/$target_sam","$output_folder/$all_deletion_outputsam")
    extract_reads_with_deletions("$output_folder/$target_sam","$output_folder/$variantpos_deletion_outputsam",[pos for (pos, _, _) in values(mut_pos_dict)])

end

function makeCIGARarray(l_CIGAR)
    l_CIGAR_num_array = split(l_CIGAR,r"[0-9]{1,3}")[2:end]
    l_CIGAR_code_array = split(l_CIGAR,r"[A-Z]")[1:end-1]
    return (l_CIGAR_num_array,l_CIGAR_code_array)
 end
 
 
 function makeCIGARstring(cigar)
 
    CIGAR = makeCIGARarray(cigar)
    CIGARstring = "" 
    for (operation, num) in zip(CIGAR[1], CIGAR[2])
 
        num = parse(Int, num)
 
        if operation =="M" || operation =="D" || operation =="I"
            CIGARstring = CIGARstring*repeat(operation,num)
        else
            CIGARstring = CIGARstring*repeat("S",num)
        end
    end
    return CIGARstring
 end

function count_D(string::AbstractString)
    count = 0
    for char in string
        if char == 'D'
            count += 1
        end
    end
    return count
end
function count_I(string::AbstractString)
    count = 0
    for char in string
        if char == 'I'
            count += 1
        end
    end
    return count
end
function count_M(string::AbstractString)
    count = 0
    for char in string
        if char != 'D' && char != 'I' 
            count += 1
        end
    end
    return count
end
function count_H(string::AbstractString)
    count = 0
    for char in string
        if char == 'H'
            count += 1
        end
    end
    return count
end
 
 function calc_mean_quality(QUAL::String)
     QUAL_LIST = Float64[]
     for q in QUAL
         push!(QUAL_LIST, Float64(q)-33.0)
     end
     return mean(QUAL_LIST)
 end
 
 
 function make_pairwise_aln(read::String, ref::String, MappedPosition::Int64, cigar::String, quality::String)

     CIGARstring = makeCIGARstring(cigar)
     pairwise_aln = [['X', 'S', 'X', '!'] for i in 1:MappedPosition-1]::Array{Array{Char,1},1}  
     counter= 1
     deletion_counter = 0
     for i in 1:length(CIGARstring) 
         if CIGARstring[i] == 'M'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, [read[j], CIGARstring[i], ref[counter+MappedPosition-1],quality[j]]) 
             counter += 1
         elseif CIGARstring[i] == 'I'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, [read[j], CIGARstring[i], '-', quality[j]])
         elseif CIGARstring[i] == 'D'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, ['-', CIGARstring[i], ref[counter+MappedPosition-1], quality[j]]) 
             counter += 1
             deletion_counter +=1
         elseif CIGARstring[i] == 'S' && counter!=1
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, ['X', CIGARstring[i], 'X', quality[j]])
         else
         end
     end
     return pairwise_aln
 end
 
 
 function create_base_call_table(ref::String)
     fulllen = length(ref)
     base_call_table = convert(DataFrame, zeros(fulllen,7))
     names!(base_call_table, [:A,:T,:G,:C,:N,:D,:I])
     return base_call_table
 end
 
 
 
 
 
 function make_single_base_call_table(ref_fasta_dir, samfile, ID;modeDI=true, qualityvalue4filter::Int64=20, meanQfilter::Float64= 20.0)
     RNAid = ""
     refseq = ""
 
     fasta = open(FASTA.Reader, ref_fasta_dir)
     for record in fasta
         RNAid = FASTA.identifier(record)
         if occursin(string(ID), RNAid)
             refseq = string(FASTA.sequence(record))
             break
         end
 
     end
     close(fasta)
     println(RNAid)
     println(refseq)
     bct = create_base_call_table(refseq)
     println("start loading ", samfile)
     f= open(samfile, "r")
     for line in eachline(f)
         if !startswith(line, "@")
             RNAID = string(split(line, "\t")[3])
             if RNAID == RNAid
                 bct = update_basecall_table(bct, refseq, line; modeDI =modeDI, qualityvalue4filter=qualityvalue4filter, meanQfilter= meanQfilter)
             end
         end
     end
     return bct
 end
 function update_basecall_table(base_call_table::DataFrame, REF::String, line::String; modeDI=modeDI::Bool, qualityvalue4filter::Int64=20, meanQfilter::Float64= 20.0, exclude_AUXMAP=true)
     nuc2numDict = Dict('A'=>1,'T'=>2,'U'=>2, 'G'=>3,'C'=>4, 'N'=>5, 'D'=>6,'I'=>7)
     FLAG_list=[]
     tokens    = split(line, "\t")
     FLAG      = parse(Int64, tokens[2])
     MAPPOS    = parse(Int64, tokens[4])
     CIGAR     = string(tokens[6])
     SEQ       = string(tokens[10])
     QUAL      = string(tokens[11])
     AUX_MAP   = string(tokens[end])
     MEAN_QUAL = calc_mean_quality(QUAL) 
     
     if (occursin("Library", AUX_MAP) && exclude_AUXMAP) || (MEAN_QUAL < meanQfilter)
         return base_call_table
     end
     push!(FLAG_list,FLAG)
 
     if !modeDI
         if occursin("D", CIGAR) || occursin("I", CIGAR)
             return base_call_table
         end
     end
 
     if occursin("H", CIGAR)
         return base_call_table
     end
 
     if FLAG == 16 || FLAG == 99 || FLAG == 147 || FLAG == 83 || FLAG == 167 || FLAG == 73 || FLAG ==  97 || FLAG == 145 || FLAG == 0
 
         pairwisealn = make_pairwise_aln(SEQ, REF, MAPPOS, CIGAR, QUAL) 
         counter = 0
         insert_counter = 0
 
         for i in 1:length(pairwisealn)
             readnucleotide, CIGARchar, elementtype, nucleotide_qc = pairwisealn[i]
             if i+1 <= length(pairwisealn)
                 nextCIGARchar = pairwisealn[i+1][2]
             else
                 nextCIGARchar = '@' 
             end
 
             if CIGARchar != 'S' 
                 counter+=1
             end
 
             if CIGARchar == 'D' && nextCIGARchar != 'D' 
                 indexnum = nuc2numDict['D']
                 if nextCIGARchar != 'D'
                     try
                         base_call_table[counter - insert_counter+MAPPOS-1, indexnum] += 1
                     catch e
                         println(e)
                         println("Error1",(readnucleotide, CIGARchar, elementtype, nucleotide_qc))
                         println(line)
                         println(pairwisealn)
                     end
                 end
 
             #2) insertion
             elseif CIGARchar == 'I' && counter - insert_counter-1>0 
 
                 if convert(Int64, nucleotide_qc)-33 >= qualityvalue4filter
                     indexnum = nuc2numDict['I']
                     if nextCIGARchar != 'I'
                         try
                             base_call_table[counter - insert_counter+MAPPOS-1, indexnum] +=1
                         catch e
                             print("Error2",(readnucleotide, CIGARchar, elementtype, nucleotide_qc))
                         end
                     end
                 end
                 insert_counter += 1
 
             elseif CIGARchar =='S' 
                 continue
 

             elseif CIGARchar == 'M' 
                 if convert(Int64, nucleotide_qc)-33 >= qualityvalue4filter
                     indexnum = nuc2numDict[readnucleotide]
                     try
                         base_call_table[counter - insert_counter+MAPPOS-1, indexnum] +=1
                     catch e

                     end
                 end
 

             end 
         end 
     end 
     
     return base_call_table
 end
 
function make_base_call_table(ref_fasta_dir, samfile, bc_output_dir; modeDI=true, qualityvalue4filter::Int64=20, meanQfilter::Float64= 20.0)

     samID = split(split(samfile,"/")[end], ".divided.sam")[1] 
     id2ref = Dict{String, String}()
     id2table = Dict{String, DataFrame}() 
     FLAG_list=[]
     fasta = open(FASTA.Reader, ref_fasta_dir)
     
     for record in fasta
         id2ref[FASTA.identifier(record)] = FASTA.sequence(record)
         id2table[FASTA.identifier(record)] = create_base_call_table(string(FASTA.sequence(record)))
     end
     close(fasta)
 
     f= open(samfile, "r")
     for line in eachline(f)
         if !startswith(line, "@")
             RNAID = string(split(line, "\t")[3])

             FLAG=string(split(line, "\t")[2])
             if haskey(id2ref, RNAID)
                 refseq = id2ref[RNAID]
                 id2table[RNAID] = update_basecall_table(id2table[RNAID], refseq, line; modeDI =modeDI, qualityvalue4filter=qualityvalue4filter, meanQfilter= meanQfilter)
             end
         end
     end
   
     close(f)
     for (RNAID, base_call_table) in id2table
         fname = "$(string(samID)).csv"
         if occursin("Ref",RNAID)
            CSV.write("$bc_output_dir/$fname",base_call_table)
         end

     end
    
 end



function make_mutationrate(call_table)
    mutationrate_array =
    map(eachrow(call_table)) do col
        (sum(col)-maximum(col))./sum(col)
    end
    replace!(mutationrate_array,NaN=>0.0)
end


function num2seq(dict::Dict, num)
    num = string(num)
    name = ""::String
    ref = ""::String

    for (k::String, v::String) in dict
        if occursin(num, k)
            name = k
            ref = v
        end
    end
    return (name, ref)::Tuple
end

function load_fasta(fasta_dir)
    dict = Dict()
    open(fasta_dir, "r") do reader
        identifier =""
        sequence=""
        for line in eachline(reader)
            if startswith(line, ">")
                identifier = replace(split(line, " ")[1], ">"=>"")
            else
                sequence =line
                dict[identifier] = sequence
            end
        end
    end
    return dict
end


function load_call_table(sampleID, seqID, suffix = "")
    fname = string(sampleID)*"_"*string(seqID)*suffix*".csv"
    df = CSV.read(fname)
    return df
end



function save_log(args_dict, log_dir)
    """Create logging directory structure according to log_dir."""
    mkpath(log_dir)
    fname = "log_"*string(now())*".txt"
    fname = replace(fname, ":"=>"-")

    open(join([log_dir, fname], "/"), "w") do f
        for (k,v) in args_dict
            write(f, string(k, ":", v, "\n"))
        end
    end
    println("Saved ", fname)
    return log_dir
end

function makeBC(fastafile,output_folder)

    samfile_directory="$output_folder"
    samfiles = [samfile for samfile in readdir("$samfile_directory") if occursin("divided.sam", samfile)]
    mkpath("$output_folder/basecalltable")
    pmap(
        samfile->make_base_call_table(
            "$fastafile",
            "$output_folder/$samfile",
            "$output_folder/basecalltable"
    ; modeDI=true, qualityvalue4filter=30, meanQfilter= 20.0),
        samfiles)
end

function main()
    args = parse_args() 
    fasta_path = args["fasta_path"]
    sam_dir=args["sam_dir"]
    output_folder=mkpath("./Demo/Output_file")
    samfiles = [samfile for samfile in readdir("$sam_dir") if occursin(".sam", samfile)]
    foreach(samfile -> begin
    make_divided_sam(fasta_path, "$sam_dir/$samfile", output_folder)
    makeBC(fasta_path, output_folder)
    end, samfiles)
end
main()

##julia ./BIVID.jl --fasta_path ./Demo/Input_file/test_FASTA_G4I8.txt --sam_dir ./Demo/Input_file/input_sam 
