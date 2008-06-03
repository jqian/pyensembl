import sys

# METHODS FOR AUTOMATIC DOWNLOADING OF RESOURCES

def copy_to_file(f,ifile=None,newpath=None,blocksize=8192000):
    'copy from file obj f to ifile (or create newpath if given)'
    if newpath is not None:
        ifile = file(newpath,'w')
    try:
        while True:
            s = f.read(blocksize)
            if s == '': break
            ifile.write(s)
    finally:
        if newpath is not None:
            ifile.close()
        f.close()

def do_gunzip(filepath,newpath=None):
    'gunzip the target using Python gzip module'
    from gzip import GzipFile
    if newpath is None:
        newpath = filepath[:-3]
    f = GzipFile(filepath)
    copy_to_file(f,newpath=newpath)
    return newpath

def run_gunzip(filepath,newpath=None):
    'run gunzip program as a sub process'
    from subprocess import call # secure way to run process
    retcode = call(['gunzip',filepath])
    if retcode != 0:
        raise OSError('gunzip "%s" failed!' % filepath)
    return filepath[:-3] # DROP THE .gz SUFFIX

def run_unzip(filepath,newpath=None,singleFile=False,**kwargs):
    '''run unzip program as a sub process,
    save to single file newpath if desired.'''
    if newpath is None:
        newpath = filepath[:-4] # DROP THE .zip SUFFIX
    from subprocess import Popen,call
    if singleFile: # concatenate all files into newpath
        import os
        ifile = file(newpath,'w')
        try:
            p = Popen(['unzip', '-p',filepath], stdout=ifile)
            pid,status = os.waitpid(p.pid, 0) # wait for unzip to exit
        finally:
            ifile.close()
    else: # just unzip the package as usual
        status = call(['unzip',filepath])
    if status != 0:
        raise OSError('unzip "%s" failed!' % filepath)
    return newpath


def create_dir_if_needed(path):
    'ensure that this directory exists, by creating it if needed'
    import os
    if os.path.isdir(path):
        return # directory exists so nothing to do
    create_dir_if_needed(os.path.dirname(path)) # ensure parent exists
    os.mkdir(path) # create this directory

def create_file_with_path(basepath, filepath):
    'create file in write mode, creating parent directory(s) if needed'
    import os.path
    newpath = os.path.join(basepath, filepath)
    create_dir_if_needed(os.path.dirname(newpath))
    return file(newpath, 'w')
        

def do_unzip(filepath, newpath=None,singleFile=False,**kwargs):
    'extract zip archive, to single file given by newpath if desired'
    # WARNING: zipfile module reads entire file into memory!
    if newpath is None:
        newpath = filepath[:-4]
    from zipfile import ZipFile
    t = ZipFile(filepath, 'r')
    try:
        if singleFile: # extract to a single file
            ifile = file(newpath,'w')
            try:
                for name in t.namelist():
                    ifile.write(t.read(name)) # may run out of memory!!
            finally:
                ifile.close()
        else: # extract a bunch of files as usual
            for name in t.namelist():
                ifile = create_file_with_path(newpath, name)
                ifile.write(t.read(name)) # may run out of memory!!
                ifile.close()
    finally:
        t.close()
    return newpath


def do_untar(filepath,mode='r|',newpath=None,singleFile=False,**kwargs):
    'extract tar archive, to single file given by newpath if desired'
    if newpath is None:
        newpath = filepath + '.out'
    import tarfile
    t = tarfile.open(filepath, mode)
    try:
        if singleFile: # extract to a single file
            ifile = file(newpath,'w')
            try:
                for name in t.getnames():
                    f = t.extractfile(name)
                    copy_to_file(f,ifile)
            finally:
                ifile.close()
        else: # extract a bunch of files as usual
            import os
            t.extractall(os.path.dirname(newpath))
    finally:
        t.close()
    return newpath


def uncompress_file(filepath,**kwargs):
    '''stub for applying appropriate uncompression based on file suffix
    (.tar .tar.gz .tgz .tar.bz2 .gz and .zip for now)'''
    if filepath.endswith('.zip'):
        print >>sys.stderr, 'unzipping %s...' % filepath
        try:
            return run_unzip(filepath,**kwargs)
        except OSError:
            return do_unzip(filepath, **kwargs)
    elif filepath.endswith('.tar'):
        print >>sys.stderr, 'untarring %s...' % filepath
        return do_untar(filepath,newpath=filepath[:-4],**kwargs)
    elif filepath.endswith('.tgz'):
        print >>sys.stderr, 'untarring %s...' % filepath
        return do_untar(filepath,mode='r:gz',newpath=filepath[:-4],**kwargs)
    elif filepath.endswith('.tar.gz'):
        print >>sys.stderr, 'untarring %s...' % filepath
        return do_untar(filepath,mode='r:gz',newpath=filepath[:-7],**kwargs)
    elif filepath.endswith('.tar.bz2'):
        print >>sys.stderr, 'untarring %s...' % filepath
        return do_untar(filepath,mode='r:bz2',newpath=filepath[:-8],**kwargs)
    elif filepath.endswith('.gz'):
        print >>sys.stderr, 'gunzipping %s...' % filepath
        try:  # could use gzip module, but it's two times slower!!
            return run_gunzip(filepath) # run as sub process
        except OSError: # on Windows, have to run as python module
            return do_gunzip(filepath)
    
    return filepath # DEFAULT: NOT COMPRESSED, SO JUST HAND BACK FILENAME

def download_monitor(bcount,bsize,totalsize):
    'show current download progress'
    bytes = bcount*bsize
    print >>sys.stderr,'downloaded %s bytes (%2.1f%%)...' \
          % (bytes,bytes*100./totalsize)

def download_unpickler(path,filename,kwargs):
    'try to download the desired file, and uncompress it if need be'
    import urllib,classutil,os
    if filename is None:
        filename = os.path.basename(path)
    filepath = os.path.join(classutil.get_env_or_cwd('PYGRDATADOWNLOAD'),filename)
    print >>sys.stderr,'Beginning download of %s to %s...' %(path,filepath)
    t = urllib.urlretrieve(path,filepath,download_monitor)
    print >>sys.stderr,'Download done.'
    filepath = uncompress_file(filepath, **kwargs) # UNCOMPRESS IF NEEDED
    o = classutil.SourceFileName(filepath) # PATH TO WHERE THIS FILE IS NOW STORED
    o._saveLocalBuild = True # MARK THIS FOR SAVING IN LOCAL PYGR.DATA
    return o
download_unpickler.__safe_for_unpickling__ = 1

class SourceURL(object):
    '''unpickling this object will trigger downloading of the desired path,
    which will be cached to PYGRDATADOWNLOAD directory if any.  The value returned
    from unpickling will simply be the path to the downloaded file, as a SourceFileName'''
    def __init__(self,path,filename=None,**kwargs):
        self.path = path
        self.kwargs = kwargs
        self.filename = filename
        if path.startswith('http:'): # make sure we can read this URL
            import httplib
            conn = httplib.HTTPConnection(path.split('/')[2])
            try:
                conn.request('GET', '/'.join([''] + path.split('/')[3:]))
                r1 = conn.getresponse()
                if r1.status!=200:
                    raise OSError('http GET failed: %d %s, %s'
                                  %(r1.status,r1.reason,path))
            finally:
                conn.close()
    def __reduce__(self):
        return (download_unpickler,(self.path,self.filename,self.kwargs))

def generic_build_unpickler(cname, args, kwargs):
    'does nothing but construct the specified klass with the specified args'
    if cname=='BlastDB':
        from seqdb import BlastDB as klass
    else:
        raise ValueError('''class name not registered for unpickling security.
Add it to pygr.downloader.generic_build_unpickler if needed: ''' + cname)
    o = klass(*args, **kwargs)
    o._saveLocalBuild = True # MARK FOR LOCAL PYGR.DATA SAVE
    return o
generic_build_unpickler.__safe_for_unpickling__ = 1

class GenericBuilder(object):
    'proxy for constructing the desired klass on unpickling'
    def __init__(self, cname, *args, **kwargs):
        self.cname = cname
        self.args = args
        self.kwargs = kwargs
    def __reduce__(self):
        return (generic_build_unpickler,(self.cname,self.args,self.kwargs))
