from predict_secretome import core
from predict_secretome import utils

if __name__ == '__main__':

    parser = utils.get_parser()
    args = parser.parse_args()

    prediction = core.predictSecretome(args)

    prediction.run_predict()

