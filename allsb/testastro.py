

import logging
import json
from email.generator import BytesGenerator
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication
from email.encoders import encode_noop
from io import BytesIO
import os.path
import time

import astropy.io.fits as fits

import requests
import http.client
import allsb.utils as U

_logger = logging.getLogger(__name__)


AN_API_KEY = 'gmbziuhfdkfnpysg'
SESSION_KEY = "pc0noq04x1gdzmc3brycrtr108vytr35"

url_base = "http://nova.astrometry.net/api"
urlx_base = "http://nova.astrometry.net"


def wcs_calibrate_astrometry_net(datafile):
    _logger.debug('shape is {shape}'.format(**datafile))
    _logger.debug('calibrate in astrometry.net')
    #astro_session()

    session = requests.Session()
    # astro_session_key(session)
    filename = datafile['filename']

    response = astro_upload_file(session, SESSION_KEY, datafile)

    url = (url_base + "/submissions/{subid}").format(**response.json())
    url_jobs = url_base + "/jobs/{}"
    url_results = url_base + "/jobs/{}/calibration/"
    url_newfits = urlx_base + "/new_fits_file/{}"
    url_axyfile = urlx_base + "/axy_file/{}"
    url_corrfile = urlx_base + "/corr_file/{}"

    while True:
        try:
            response = requests.get(url)
            mm = response.json()
            print(mm)
            time.sleep(5)
            print('jobs:', mm['jobs'])
            if mm['jobs']:
                for jid in mm['jobs']:
                    if jid is not None:
                        response_job = requests.get(url_jobs.format(jid))
                        print('job', jid, response_job.text)
            print('job_calibrations:', mm['job_calibrations'])
            if mm['job_calibrations']:
                break
        except http.client.RemoteDisconnected as error:
            print(error)

    response = requests.get(url_results.format(jid))
    nfile = U.insert_adj_filename(filename, 'results', ext='.json')
    _logger.debug('save astrometry.net results in %s', nfile)
    with open(nfile, 'wb') as fd:
        fd.write(response.content)


    urls = [url_newfits.format(jid), url_corrfile.format(jid)]
    ajs_s = ['newfits', 'corrfile']

    for url, adj in zip(urls, ajs_s):
        response = requests.get(url)
        nfile = U.insert_adj_filename(filename, adj)
        _logger.debug('save in %s', nfile)
        with open(nfile, 'wb') as fd:
            fd.write(response.content)


def astro_session_key(session):

    rjson = json.dumps({"apikey": AN_API_KEY})
    result = session.post('http://nova.astrometry.net/api/login', data={'request-json': rjson})
    print(result)
    return


def cut_center(filename, center, hsize=500):
    with fits.open(filename) as hdulist:
        data = hdulist[0].data
        newdata = data[center[0] - hsize:center[0] + hsize,
                  center[1] - hsize:center[1] + hsize]
        newhdu = fits.PrimaryHDU(newdata, hdulist[0].header)
        return newhdu


def calc_filename_center(filename):
    return U.insert_adj_filename(filename, 'center')


def astro_upload_file(session, session_key, datafile):
    upload_args = {
        'scale_upper': 90.0, 'publicly_visible': 'y',
        'allow_modifications': 'd', 'scale_type': 'ul',
        'allow_commercial_use': 'd', 'scale_lower': 25.0
    }
    upload_args['session'] = session_key

    filename = datafile['filename']
    # cut filename
    res = datafile['res']

    center = (int(res[1]), int(res[0]))
    _logger.debug('computing cut image around x=%s, y=%s', center[1], center[0])
    cut_hdu = cut_center(filename, center)
    cut_filename = calc_filename_center(filename)
    cut_hdu.writeto(cut_filename, overwrite=True)

    try:
        with open(cut_filename, 'rb') as f:
            file_args = (cut_filename, f.read())
    except IOError as error:
        print(error)
        raise

    headers, data = data_encode(upload_args, file_args)
    response = session.post('http://nova.astrometry.net/api/upload', headers=headers, data=data)
    return response


class MyBytesGenerator(BytesGenerator):
    def __init__(self, fp, root=True):
        super().__init__(fp, mangle_from_=False,maxheaderlen=0)
        self.root = root

    def _write_headers(self, msg):
        # We don't want to write the top-level headers;
        # they go into Request(headers) instead.
        if self.root:
            return
        # We need to use \r\n line-terminator, but Generator
        # doesn't provide the flexibility to override, so we
        # have to copy-n-paste-n-modify.
        for h, v in msg.items():
            self._fp.write(('%s: %s\r\n' % (h,v)).encode())
        # A blank line always separates headers from body
        self._fp.write('\r\n'.encode())

    # The _write_multipart method calls "clone" for the
    # subparts.  We hijack that, setting root=False
    def clone(self, fp):
        return MyBytesGenerator(fp, root=False)


def data_encode(args, file_args):

    # Make a custom generator to format it the way we need.
    json_args = json.dumps(args)
    # If we're sending a file, format a multipart/form-data

    m1 = MIMEBase('text', 'plain')
    m1.add_header('Content-disposition',
                  'form-data; name="request-json"')
    m1.set_payload(json_args)
    m2 = MIMEApplication(file_args[1], 'octet-stream', encode_noop)
    m2.add_header('Content-disposition',
                  'form-data; name="file"; filename="%s"' % file_args[0])
    mp = MIMEMultipart('form-data', None, [m1, m2])

    fp = BytesIO()
    g = MyBytesGenerator(fp)
    g.flatten(mp)
    data = fp.getvalue()
    headers = {'Content-type': mp.get('Content-type')}

    return headers, data
