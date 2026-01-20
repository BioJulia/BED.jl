module FileIOExt

using BED, FileIO

function BED.fileio_load(filepath::File{format"BED"})
    
    reader = open(BED.Reader, FileIO.filename(filepath))
    records = try 
        records = [record for record in reader]
    catch err
        close(reader)
        rethrow(err)
    end
    close(reader)

    records
end

function BED.fileio_save(filepath::File{format"BED"}, records)

    writer = open(BED.Writer, FileIO.filename(filepath))
    try 
        for record in records
            write(writer, record)
        end
    catch err
        close(writer)
        rethrow(err)
    end
    close(writer)
    filepath
end

end