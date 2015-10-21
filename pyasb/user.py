import argparse

from pyasb.info import __version__

config_file_default = 'config.cfg'

def main(args=None):

    thisversion = '%(prog)s {}'.format(__version__)

    parser = argparse.ArgumentParser(prog='PyASB',
                                     description='Welcome to %(prog)s (Python All-Sky Brightness pipeline)'
                                     )
    parser.add_argument('--version', action='version', version=thisversion)
    parser.add_argument('-i', metavar="input_allsky_image", nargs='+',
                        dest='fits_filename_list',
                        help='All Sky image you want to be analyzed')
    parser.add_argument('-c', action='store', dest='configfile', metavar="config_file",
                        default=config_file_default,
                        help='Use alternative config file')
    parser.add_argument('-d', metavar="[year,month,day]",
                        help='Date to be analyzed (AstMon-UCM only), '
                             'month and day are optionaUse alternative config file')
    parser.add_argument('-om', action='store_true',
                        help='Output star map image, full or relative path, '
                             'if no output file, show the map on screen')
    parser.add_argument('-ot', metavar="output_photometric_table path",
                        help='Output photometric table full or relative path, '
                             'if no output file, show the table on screen')
    parser.add_argument('-or', metavar="output_results_summary path",
                        help='Summary of analysis, fit parameters and zenith SB '
                             'full or relative path. If no output file, '
                             'show the table on screen')
    parser.add_argument('-ob', metavar="output_bouguerfit_graph path",
                        help='Output bouguer-fit graph, full or relative path. '
                             'If no output file, show the graph on screen')
    parser.add_argument('-ocm', metavar="output_cloudmap_image path",
                        help='Output cloud map image, full or relative path, '
                             'if no output file, show the map on screen')
    parser.add_argument('-oct', metavar="output_cloudmap_image path",
                        help='Output cloud map image, full or relative path, '
                             'if no output file, show the map on screen')
    parser.add_argument('-os', metavar="output_skybrightness_graph path",
                        help='Output Sky Brightness graph, full or relative path,'
                             ' if no output file, show the graph on screen')
    parser.add_argument('-ost', metavar="output_skybrightness_table path",
                        help='Output Sky Brightness table, full or relative path, '
                             'if no output file, show the graph on screen')
    parser.add_argument('--path',
                        help='save path')
    resargs = parser.parse_args(args)
    return resargs


if __name__ == '__main__':
    main()
