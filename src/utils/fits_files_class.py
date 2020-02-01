

class FitsClass():
    def __init__(self, filepath, name, offset_data, size):
        self.file_path = filepath 
        self.open_file = open(filepath, mode = 'rb')

        self.name = name 

        self.offset_data = offset_data 
        self.size = size 


    def __del__(self):
        """
            Makes sure that the file is closed
        """
        self.open_file.close()

    def seek(self, seek_pos):
        self.open_file.seek(int(seek_pos))

    @property 
    def file(self):
        return self.open_file


if __name__ == '__main__':
    import tarfile 
    s = '/home/amiguel/SERVALHOME/serval/data/HARPS/gj699/ADP.2014-09-17T11-23-26.740.tar'
    tar = tarfile.open(s)

    for member in tar.getmembers():
      if '_e2ds_A.fits' in member.name: 
        extr = member

        s = FitsClass(filepath = s, 
                offset_data = extr.offset_data,
                name = extr.name,
                size=extr.size) 


    s.seek(s.offset_data)    
    args = (b'INSTRUME', b'OBJECT', b'MJD-OBS')
    for card in iter(lambda:s.file.read(80), ''):   # read in 80 byte blocks
        if card.startswith(args):
            print('yes')

