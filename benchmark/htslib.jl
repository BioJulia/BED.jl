# Minimal htslib bindings
# =======================
#
# Thin `ccall` wrappers around libhts from htslib_jll, covering just enough of the API to compress, index, scan and query BED files.
# The structs mirror their C counterparts in htslib 1.10 or later, where `hts_pos_t` is a 64-bit integer.
# Records are read with `tbx_readrec`, the reader that tabix itself uses: it reads a line and parses its sequence name, begin and end.

module HTSlib

using Libdl: dlsym
using htslib_jll: htslib_jll, libhts, libhts_handle

# kstring_t
struct KString
    length::Csize_t
    capacity::Csize_t
    data::Ptr{UInt8}
end

# tbx_conf_t
struct TabixConfiguration
    preset::Int32
    column_sequence::Int32
    column_begin::Int32
    column_end::Int32
    meta_char::Int32
    line_skip::Int32
end

# The leading fields of tbx_t.
struct Tabix
    configuration::TabixConfiguration
    index::Ptr{Cvoid}
    dictionary::Ptr{Cvoid}
end

# The leading fields of hts_itr_t, up to the interval of the record last read.
struct Iterator
    flags::UInt32
    tid::Cint
    offset_count::Cint
    i::Cint
    region_count::Cint
    start_position::Int64
    end_position::Int64
    region_list::Ptr{Cvoid}
    current_tid::Cint
    current_region::Cint
    current_interval::Cint
    current_start_position::Int64
    current_end_position::Int64
end

tabix_configuration_bed() = cglobal((:tbx_conf_bed, libhts), TabixConfiguration)

"""
    bgzip(source, destination)

Compress the file `source` into the BGZF file `destination` with htslib.
"""
function bgzip(source::AbstractString, destination::AbstractString)
    data = read(source)
    file = ccall((:bgzf_open, libhts), Ptr{Cvoid}, (Cstring, Cstring), destination, "w")
    file == C_NULL && error("bgzf_open failed for $(destination)")
    written = ccall((:bgzf_write, libhts), Cssize_t, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t), file, data, length(data))
    written == length(data) || error("bgzf_write failed for $(destination)")
    ccall((:bgzf_close, libhts), Cint, (Ptr{Cvoid},), file) == 0 || error("bgzf_close failed for $(destination)")
    return destination
end

"""
    tabix_index(path)

Build the tabix index `path * ".tbi"` for the BGZF-compressed BED file `path`.
"""
function tabix_index(path::AbstractString)
    status = ccall((:tbx_index_build, libhts), Cint, (Cstring, Cint, Ptr{TabixConfiguration}), path, 0, tabix_configuration_bed())
    status == 0 || error("tbx_index_build failed for $(path)")
    return path * ".tbi"
end

"""
    IndexedFile(path)

An open BGZF-compressed BED file together with its loaded tabix index.
Call `close` when finished.
"""
mutable struct IndexedFile
    file::Ptr{Cvoid}
    tabix::Ptr{Tabix}
end

function IndexedFile(path::AbstractString)
    file = ccall((:hts_open, libhts), Ptr{Cvoid}, (Cstring, Cstring), path, "r")
    file == C_NULL && error("hts_open failed for $(path)")
    tabix = ccall((:tbx_index_load, libhts), Ptr{Tabix}, (Cstring,), path)
    tabix == C_NULL && error("tbx_index_load failed for $(path)")
    return IndexedFile(file, tabix)
end

function Base.close(indexed::IndexedFile)
    if indexed.tabix != C_NULL
        ccall((:tbx_destroy, libhts), Cvoid, (Ptr{Tabix},), indexed.tabix)
        indexed.tabix = C_NULL
    end
    if indexed.file != C_NULL
        ccall((:hts_close, libhts), Cint, (Ptr{Cvoid},), indexed.file)
        indexed.file = C_NULL
    end
    return nothing
end

"""
    scan(path, indexed)

Read every record of the BED file `path`, which may be plain or BGZF-compressed, with `tbx_readrec`.
The sequence names are resolved against the index of `indexed`, which must describe the same records.
Returns the number of records and the sum of their interval lengths.
"""
function scan(path::AbstractString, indexed::IndexedFile)
    # bgzf_open reads uncompressed files transparently.
    file = ccall((:bgzf_open, libhts), Ptr{Cvoid}, (Cstring, Cstring), path, "r")
    file == C_NULL && error("bgzf_open failed for $(path)")
    line = Ref(KString(0, 0, C_NULL))
    tid = Ref{Cint}(0)
    start_position = Ref{Int64}(0)
    end_position = Ref{Int64}(0)
    count = 0
    total_length = 0
    while true
        status = ccall(
            (:tbx_readrec, libhts), Cint,
            (Ptr{Cvoid}, Ptr{Tabix}, Ref{KString}, Ref{Cint}, Ref{Int64}, Ref{Int64}),
            file, indexed.tabix, line, tid, start_position, end_position,
        )
        status == -1 && break
        status < 0 && error("tbx_readrec failed on record $(count + 1) of $(path)")
        count += 1
        total_length += end_position[] - start_position[]
    end
    Libc.free(line[].data)
    ccall((:bgzf_close, libhts), Cint, (Ptr{Cvoid},), file)
    return count, total_length
end

"""
    query(indexed, region)

Iterate over the records of `indexed` overlapping `region` (for example `"chr1:1001-2000"`), as `tabix file region` does.
Returns the number of records and the sum of their interval lengths.
"""
function query(indexed::IndexedFile, region::AbstractString)
    # tbx_itr_querys and tbx_itr_next are C macros, so expand them here.
    index = unsafe_load(indexed.tabix).index
    iterator = ccall(
        (:hts_itr_querys, libhts), Ptr{Iterator},
        (Ptr{Cvoid}, Cstring, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}),
        index, region, dlsym(libhts_handle, :tbx_name2id), indexed.tabix, dlsym(libhts_handle, :hts_itr_query), dlsym(libhts_handle, :tbx_readrec),
    )
    iterator == C_NULL && error("hts_itr_querys failed for $(region)")
    stream = ccall((:hts_get_bgzfp, libhts), Ptr{Cvoid}, (Ptr{Cvoid},), indexed.file)
    line = Ref(KString(0, 0, C_NULL))
    count = 0
    total_length = 0
    while true
        status = ccall((:hts_itr_next, libhts), Cint, (Ptr{Cvoid}, Ptr{Iterator}, Ref{KString}, Ptr{Tabix}), stream, iterator, line, indexed.tabix)
        status == -1 && break
        status < 0 && error("hts_itr_next failed for $(region)")
        current = unsafe_load(iterator)
        count += 1
        total_length += current.current_end_position - current.current_start_position
    end
    Libc.free(line[].data)
    ccall((:hts_itr_destroy, libhts), Cvoid, (Ptr{Iterator},), iterator)
    return count, total_length
end

end # module
