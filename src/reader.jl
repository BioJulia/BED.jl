# BED Reader
# ==========
function appendfrom!(dst, dpos, src, spos, n)
    if length(dst) < dpos + n - 1
        resize!(dst, dpos + n - 1)
    end
    unsafe_copyto!(dst, dpos, src, spos, n)
    return dst
end

mutable struct Reader <: BioGenerics.IO.AbstractReader
    state::BioGenerics.Automa.State
    index::Union{Indexes.Tabix,Nothing}

    function Reader(stream::TranscodingStream, index=nothing)
        return new(BioGenerics.Automa.State(stream, 1, 1, false), index)
    end
end

"""
    BED.Reader(input::IO; index=nothing)
    BED.Reader(input::AbstractString; index=:auto)

Create a data reader of the BED file format.

The first argument specifies the data source.
When it is a filepath that ends with *.bgz*, it is considered to be block compression file format (BGZF) and the function will try to find a tabix index file (<filename>.tbi) and read it if any.
See <http://www.htslib.org/doc/tabix.html> for bgzip and tabix tools.

Arguments
---------
- `input`: data source
- `index`: path to a tabix file
"""
function Reader(input::IO; index=nothing)
    if isa(index, AbstractString)
        index = Indexes.Tabix(index)
    end

    stream = TranscodingStreams.NoopStream(input)

    return Reader(stream, index)

end

function Reader(filepath::AbstractString; index=:auto)
    if isa(index, Symbol) && index != :auto
        throw(ArgumentError("invalid index argument: ':$(index)'"))
    end
    if endswith(filepath, ".bgz")
        input = BGZFStreams.BGZFStream(filepath)
        if index == :auto
            index = Indexes.findtabix(filepath)
        end
    else
        input = open(filepath)
    end
    return Reader(input, index = index)
end

function Base.eltype(::Type{Reader})
    return Record
end

function BioGenerics.IO.stream(reader::Reader)
    return reader.state.stream
end

function Base.iterate(reader::Reader, nextone = Record())
    if BioGenerics.IO.tryread!(reader, nextone) === nothing
        return nothing
    end
    return copy(nextone), empty!(nextone) # Empty record for inplace reading and reuse of array allocations.
end

function GenomicFeatures.eachoverlap(reader::Reader, interval::GenomicFeatures.Interval)
    if reader.index === nothing
        throw(ArgumentError("index is null"))
    end
    return Indexes.TabixOverlapIterator(reader, interval)
end

const record_machine, file_machine = let
    alt = Automa.RegExp.alt
    cat = Automa.RegExp.cat
    rep = Automa.RegExp.rep
    opt = Automa.RegExp.opt

    record = let
        # Note: The name cannot end with a #, nor whitespace, because if so,
        # it is undistinguishable from whitespace and comments
        chrom = onexit!(onenter!(re"" | re"[^# \t\v\r\n\f][ -~]*", :pos), :record_chrom)
        chromstart = onexit!(onenter!(re"[0-9]+", :pos), :record_chromstart)
        chromend = onexit!(onenter!(re"[0-9]+", :pos), :record_chromend)
        name = onexit!(onenter!(re"[ -~]*", :pos), :record_name)
        score = onexit!(onenter!(re"[0-9]+", :pos), :record_score)
        strand = onenter!(re"[+\-.?]", :record_strand) #Note: single byte.
        thickstart = onexit!(onenter!(re"[0-9]+", :pos), :record_thickstart)
        thickend = onexit!(onenter!(re"[0-9]+", :pos), :record_thickend)
        itemrgb = cat(re"[0-9]+", opt(cat(',', re"[0-9]+", ',', re"[0-9]+")))
        onexit!(onenter!(itemrgb, :pos), :record_itemrgb)
        blockcount = onexit!(onenter!(re"[0-9]+", :pos), :record_blockcount)

        # comma-separated values
        csv(x) = cat(rep(cat(x, ',')), opt(x))

        blocksizes = csv(onexit!(onenter!(re"[0-9]+", :pos), :record_blocksizes_blocksize))
        onexit!(blocksizes, :record_blocksizes)
        blockstarts = csv(onexit!(onenter!(re"[0-9]+", :pos), :record_blockstarts_blockstart))
        onexit!(blockstarts, :record_blockstarts)

        cat(
            chrom, '\t',
            chromstart, '\t',
            chromend,
            opt(cat('\t', name,
            opt(cat('\t', score,
            opt(cat('\t', strand,
            opt(cat('\t', thickstart,
            opt(cat('\t', thickend,
            opt(cat('\t', itemrgb,
            opt(cat('\t', blockcount,
            opt(cat('\t', blocksizes,
            opt(cat('\t', blockstarts)))))))))))))))))))
    end
    onexit!(onenter!(record, :mark), :record)
    hspace = re"[ \t\v\f]"

    # Note: We cannot start with a tab, because then it is indistinguishable
    # from an empty chrom field
    blankline = opt(re"[ \v\f]" * rep(hspace))
    comment = re"#[^\n]*"

    newline = let
        lf = onenter!(re"\n", :countline)
        cat(opt('\r'), lf)
    end

    file = rep(alt(
        cat(record, newline),
        cat(blankline, newline),
        cat(comment, newline),
    ))

    map(Automa.compile, (record, file))
end

const record_actions = Dict(
    :mark => :(@mark),
    :pos => :(pos = @relpos(p)),
    :record_chrom => :(record.chrom = (pos:@relpos(p-1)); record.ncols += 1),
    :record_chromstart => :(record.chromstart = (pos:@relpos(p-1)); record.ncols += 1),
    :record_chromend => :(record.chromend = (pos:@relpos(p-1)); record.ncols += 1),
    :record_name => :(record.name = (pos:@relpos(p-1)); record.ncols += 1),
    :record_score => :(record.score = (pos:@relpos(p-1)); record.ncols += 1),
    :record_strand => :(record.strand = @relpos(p); record.ncols += 1),
    :record_thickstart => :(record.thickstart = (pos:@relpos(p-1)); record.ncols += 1),
    :record_thickend => :(record.thickend = (pos:@relpos(p-1)); record.ncols += 1),
    :record_itemrgb => :(record.itemrgb = (pos:@relpos(p-1)); record.ncols += 1),
    :record_blockcount => :(record.blockcount = (pos:@relpos(p-1)); record.ncols += 1),
    :record_blocksizes_blocksize => :(push!(record.blocksizes, (pos:@relpos(p-1)))),
    :record_blocksizes => :(record.ncols += 1),
    :record_blockstarts_blockstart => :(push!(record.blockstarts, (pos:@relpos(p-1)))),
    :record_blockstarts => :(record.ncols += 1),
    :record => :(record.filled = 1:@relpos(p-1))
)

Automa.generate_reader(
    :index!,
    record_machine,
    arguments = (:(record::Record),),
    actions = record_actions,
    # context = :(),
    initcode = :(pos = 0),
    # loopcode = :()
    # returncode = :()
) |> eval


const initcode = quote
    pos = 0
    linenum = 0
    found_record=false
    # empty!(record)
    cs, linenum = state
end

const loopcode = quote
    if found_record
        @goto __return__
    end
end

Automa.generate_reader(
    :readrecord!,
    file_machine,
    arguments = (:(record::Record), :(state::Tuple{Int,Int})),
    actions = merge(record_actions, Dict(
        :record => quote
            appendfrom!(record.data, 1, data, @markpos, p-@markpos)
            record.filled = 1:(p-@markpos)
            found_record = true
            @escape
        end,
        :countline => :(linenum += 1),
    )),
    initcode = initcode,
    loopcode = loopcode,
    returncode = :(return cs, linenum, found_record)
) |> eval


function index!(record::Record)
    stream = TranscodingStreams.NoopStream(IOBuffer(record.data))
    cs = index!(stream, record)
    if cs != 0
        throw(ArgumentError("Invalid BED record. Machine failed to transition from state $(cs)."))
    end
    return record
end

"""
    read!(rdr::Reader, rec::Record)
Read a `Record` into `rec`; overwriting or adding to existing field values.
It is assumed that `rec` is already initialized or empty.
"""
function Base.read!(rdr::Reader, record::Record)

    cs, ln, found = readrecord!(rdr.state.stream, record, (rdr.state.state, rdr.state.linenum))

    rdr.state.state = cs
    rdr.state.linenum = ln
    rdr.state.filled = found

    if found
        return record
    end

    if cs == 0 || eof(rdr.state.stream)
        throw(EOFError())
    end

    throw(ArgumentError("Malformed BED file record at line $(ln). Machine failed to transition from state $(cs)."))
end
