
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
         "--output_dir"
            help = "Output dir for divided samfiles"
            default = "./Demo/Output_file/"
    end
    return ArgParse.parse_args(s)  
end
 

struct VariantGroup
    ref_id::String
    ref_num::String
    ref_seq::String
    gene::String
    variants::Vector{Tuple{String,Int,Char,Char}}
end

struct VariantDef
    num::String
    pos::Int
    minor::Char
    major::Char
    path::String
end

function read_fasta_dict(path::AbstractString)
    entries = Dict{String,String}()
    open(path, "r") do io
        header = nothing
        seq_buf = IOBuffer()
        for raw in eachline(io)
            line = strip(raw)
            isempty(line) && continue
            if startswith(line, ">")
                if header !== nothing
                    entries[header] = String(take!(seq_buf))
                end
                header = strip(line[2:end])
            else
                write(seq_buf, line)
            end
        end
        if header !== nothing
            entries[header] = String(take!(seq_buf))
        end
    end
    return entries
end

parse_header(header::String) = split(header, '_', limit=3)

function build_variant_groups(path::AbstractString)
    entries = read_fasta_dict(path)
    references = Dict{String,Tuple{String,String}}()
    for (header, seq) in entries
        parts = parse_header(header)
        length(parts) >= 2 || error("Unexpected FASTA header format: $header")
        gene = parts[2]
        if occursin("Reference", header)
            references[gene] = (header, seq)
        end
    end
    isempty(references) && error("No reference entries found in FASTA")

    groups = VariantGroup[]
    for (gene, (ref_id, ref_seq)) in references
        ref_parts = parse_header(ref_id)
        ref_num = ref_parts[1]
        variants = Tuple{String,Int,Char,Char}[]
        for (header, seq) in entries
            occursin("Reference", header) && continue
            parts = parse_header(header)
            length(parts) >= 2 || continue
            parts[2] == gene || continue
            length(seq) == length(ref_seq) || error("Length mismatch between $header and $ref_id")
            diffs = [(i, ref_seq[i], seq[i]) for i in eachindex(seq) if ref_seq[i] != seq[i]]
            isempty(diffs) && continue
            length(diffs) == 1 || error("Multiple differences detected for $header; expected single SNP")
            pos, maj, min = diffs[1]
            push!(variants, (parts[1], pos, uppercase(min), uppercase(maj)))
        end
        push!(groups, VariantGroup(ref_id, ref_num, ref_seq, gene, variants))
    end
    return groups
end

function read_headers_and_first_alignment(io)
    headers = String[]
    first_alignment::Union{Nothing,String} = nothing
    while !eof(io)
        line = readline(io)
        if startswith(line, '@')
            push!(headers, line)
        else
            first_alignment = line
            break
        end
    end
    return headers, first_alignment
end

function write_headers!(io::IO, headers::Vector{String})
    for header in headers
        write(io, header)
        endswith(header, "\n") || write(io, '\n')
    end
end

function ensure_writer!(writers::Dict{String,IO}, path::String, headers::Vector{String})
    if !haskey(writers, path)
        isfile(path) && rm(path)
        io = open(path, "w")
        write_headers!(io, headers)
        writers[path] = io
    end
    return writers[path]
end

function query_base_at_reference(fields::Vector{SubString{String}}, seq_chars::Vector{Char}, target_pos::Int)
    cigar = fields[6]
    cigar == "*" && return nothing
    ref_start = parse(Int, fields[4])
    ref_pos = ref_start
    read_idx = 1
    for match in eachmatch(r"(\d+)([MIDNSHP=X])", cigar)
        len = parse(Int, match.captures[1])
        op = match.captures[2]
        if op == "M" || op == "=" || op == "X"
            for _ in 1:len
                if ref_pos == target_pos
                    return read_idx <= length(seq_chars) ? seq_chars[read_idx] : nothing
                end
                ref_pos += 1
                read_idx += 1
            end
        elseif op == "I" || op == "S"
            read_idx += len
        elseif op == "D" || op == "N"
            if target_pos >= ref_pos && target_pos < ref_pos + len
                return nothing
            end
            ref_pos += len
        end
    end
    return nothing
end

function classify_alignment(fields::Vector{SubString{String}}, seq_chars::Vector{Char}, variant_defs::Vector{VariantDef})
    matched_variant::Union{Nothing,VariantDef} = nothing
    for variant in variant_defs
        base = query_base_at_reference(fields, seq_chars, variant.pos)
        base === nothing && return nothing
        nb = uppercase(base)
        if nb == variant.minor
            matched_variant === nothing || return nothing
            matched_variant = variant
        elseif nb == variant.major
            continue
        else
            return nothing
        end
    end
    return matched_variant === nothing ? :major : matched_variant
end

function extract_reads_by_fasta_id(input_samfile::AbstractString, output_samfile::AbstractString, fasta_id::AbstractString)
    open(input_samfile, "r") do input_file
        headers, first_alignment = read_headers_and_first_alignment(input_file)
        writer = open(output_samfile, "w")
        write_headers!(writer, headers)
        function maybe_write(line)
            fields = split(line, '\t')
            if fields[3] == fasta_id
                write(writer, line)
                write(writer, '\n')
            end
        end
        if first_alignment !== nothing
            maybe_write(first_alignment)
        end
        while !eof(input_file)
            maybe_write(readline(input_file))
        end
        close(writer)
    end
end

function partition_target_sam!(temp_path::AbstractString, ref_path::AbstractString,
                               all_path::AbstractString, unassigned_path::AbstractString,
                               variant_defs::Vector{VariantDef})
    open(temp_path, "r") do input
        headers, first_alignment = read_headers_and_first_alignment(input)
        writers = Dict{String,IO}()
        get_writer(path::String) = ensure_writer!(writers, path, headers)
        all_io = get_writer(all_path)
        function dispatch(line)
            fields = split(line, "\t")
            length(fields) >= 10 || return
            seq_chars = collect(fields[10])
            result = classify_alignment(fields, seq_chars, variant_defs)
            line_out = string(line, "\n")
            if result === nothing
                unassigned_io = get_writer(unassigned_path)
                write(unassigned_io, line_out)
            else
                target_io = result === :major ? get_writer(ref_path) : get_writer(result.path)
                write(target_io, line_out)
            end
            write(all_io, line_out)
        end
        if first_alignment !== nothing
            dispatch(first_alignment)
        end
        while !eof(input)
            dispatch(readline(input))
        end
        for io in values(writers)
            close(io)
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

            if occursin("D",cigar_string)==true
                write(output_file, line, "\n")
            end
        end
    end
end
end
function make_divided_sam(fasta_path::AbstractString, sam_path::AbstractString,output_dir::AbstractString)

    groups = build_variant_groups(fasta_path)
    mkpath(output_dir)
    Parent_sam=sam_path
    sample_label = sample_label = splitext(basename(Parent_sam))[1]
    sample_all_path = joinpath(output_dir, string(sample_label, "_all.sam"))
    isfile(sample_all_path) && rm(sample_all_path)
    for group in groups
        gene = group.gene
        println(sample_label)
        println(gene)
        println(output_dir)
        ref_path = joinpath(output_dir, string(sample_label, "_", group.ref_num, "_", gene, ".sam"))
        all_path = joinpath(output_dir, string(sample_label, "_all_", gene, ".sam"))
        temp_path = joinpath(output_dir, string(sample_label, "_", gene, "_tmp.sam"))
        unassigned_path = joinpath(output_dir, string(sample_label, "_unassigned_", gene, ".sam"))
        isfile(temp_path) && rm(temp_path)
        println(ref_path)
        isfile(unassigned_path) && rm(unassigned_path)
        extract_reads_by_fasta_id(Parent_sam, temp_path, group.ref_id)
        if isempty(group.variants)
            isfile(ref_path) && rm(ref_path)
            mv(temp_path, ref_path; force=true)
            cp(ref_path, all_path; force=true)
            append_sample_all(sample_all_path, ref_path)
            continue
        end
        variant_defs = VariantDef[]
        for (num, pos, minor, major) in group.variants

            variant_path = joinpath(output_dir, string(sample_label, "_", num, "_", gene, ".sam"))
            println(variant_path)
            isfile(variant_path) && rm(variant_path)
            push!(variant_defs, VariantDef(num, pos, minor, major, variant_path))
        end
        isfile(ref_path) && rm(ref_path)
        isfile(all_path) && rm(all_path)
        partition_target_sam!(temp_path, ref_path, all_path, unassigned_path, variant_defs)
        rm(temp_path; force=true)
        append_sample_all(sample_all_path, all_path)
        extract_reads_with_alldeletions(all_path,replace(all_path,".sam"=>".deletion.sam"))
        extract_reads_with_alldeletions(unassigned_path,replace(all_path,".sam"=>".deletion.sam"))
    end

    return nothing
end

function append_sample_all(sample_all_path::AbstractString, gene_all_path::AbstractString; include_header=false)
    if !isfile(sample_all_path)
        if include_header
            cp(gene_all_path, sample_all_path; force=true)
        else
            open(gene_all_path, "r") do fin
                open(sample_all_path, "w") do fout
                    for line in eachline(fin)
                        startswith(line, "@") && continue
                        write(fout, line)
                        endswith(line, "\n") || write(fout, "\n")
                    end
                end
            end
        end
        return
    end
    open(sample_all_path, "a") do fout
        open(gene_all_path, "r") do fin
            skip_header = !include_header
            for line in eachline(fin)
                if skip_header && startswith(line, "@")
                    continue
                end
                write(fout, line)
                endswith(line, "\n") || write(fout, "\n")
            end
        end
    end
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
 
 
 function calc_mean_quality(QUAL::String)
     QUAL_LIST = Float64[]
     for q in QUAL
         push!(QUAL_LIST, Float64(q)-33.0)
     end
     return mean(QUAL_LIST)
 end
 
 
 function make_pairwise_aln(read::String, ref::String, MappedPosition::Int64, cigar::String, quality::String)

     CIGARstring = makeCIGARstring(cigar)
 

     pairwise_aln = [['X', 'S', 'X', '!'] for i in 1:MappedPosition-1]::Array{Array{Char,1},1}  #11/26に付け足した. これが無いと逆鎖が大変なことになる.
     #TODO: we should pre-allocate memory for this array.
 
     counter= 1

     deletion_counter = 0
     for i in 1:length(CIGARstring) #iはcigar stringのloc
         if CIGARstring[i] == 'M'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, [read[j], CIGARstring[i], ref[counter+MappedPosition-1],quality[j]]) #mapping posiの分だけずらす.
             counter += 1
         elseif CIGARstring[i] == 'I'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, [read[j], CIGARstring[i], '-', quality[j]])
         elseif CIGARstring[i] == 'D'
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, ['-', CIGARstring[i], ref[counter+MappedPosition-1], quality[j]])  #mapping posiの分だけずらす.
             counter += 1
             deletion_counter +=1
         elseif CIGARstring[i] == 'S' && counter!=1 #exclude the case cigar starts with S
             j = i-deletion_counter
             pairwise_aln = push!(pairwise_aln, ['X', CIGARstring[i], 'X', quality[j]])
         else
             #print(CIGARstring[i])
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
 
 
 
 println("base_call_table.jl installed.")
 
 
 
  function make_base_call_table(ref_fasta_dir, samfile, outputdir; modeDI=true, qualityvalue4filter::Int64=20, meanQfilter::Float64= 20.0)
 
     """
     https://samtools.github.io/hts-specs/SAMv1.pdf
     samfile: YYYYMMDD_N_250merC.sam (N=samID)
     in the fasta file, identifier should be "NAME:NUMBER_DESCRIPTION". Then this function pick up the information.
     """
     samID = split(split(samfile,"/")[end], ".sam")[1] 
     id2ref = Dict{String, String}()
     id2table = Dict{String, DataFrame}() #my library contains 1800seq.
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
             if haskey(id2ref, RNAID) #exclude "*"
 
                 refseq = id2ref[RNAID]
     
                 id2table[RNAID] = update_basecall_table(id2table[RNAID], refseq, line; modeDI =modeDI, qualityvalue4filter=qualityvalue4filter, meanQfilter= meanQfilter)
             end
         end
     end
   
     close(f)
     for (RNAID, base_call_table) in id2table
         
         fname = outputdir*"/"*string(samID)*".csv"
         
         CSV.write(fname,base_call_table)
     end
    
 end
 using DataFrames
using CSV
using Statistics
using Dates


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

    #番号の検索
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

function makeBC(fastafile,output_dir)

    samfile_directory="$output_dir"
    samfiles = [samfile for samfile in readdir("$samfile_directory") if occursin(".sam", samfile)]
    mkpath("$output_dir/basecalltable")
    pmap(
        samfile->make_base_call_table(
            "$fastafile",
            "$output_dir/$samfile",
            "$output_dir/basecalltable"
    ; modeDI=true, qualityvalue4filter=30, meanQfilter= 20.0),
        samfiles)
end
function main()
    args = parse_args() 
    fasta_path = args["fasta_path"]
    sam_dir=args["sam_dir"]
    output_dir=args["output_dir"]
    samfiles = [samfile for samfile in readdir("$sam_dir") if occursin(".sam", samfile)]
    foreach(samfile -> begin
    make_divided_sam(fasta_path, "$sam_dir/$samfile", output_dir)
    println(samfile)
    makeBC(fasta_path, output_dir)
    end, samfiles)
end
main()