<html>
    <body>
        <script>

            function sigmoid(x) {
                return 1 / (1 + Math.exp(-x));
            }

            function sigmoid2(x) {
                return 2 * sigmoid(x * 4.9) - 1;
            }

            var w1 = 0.3;
            var w2 = 0.2;

            var x1 = -0.71;
            var x2 = 0.321;
            d_sigmoid = function (x) {
                return sigmoid(x) * (1 - sigmoid(x));
            };
            var sum = w1 * x1 + w2 * x2;
            var y1 = sigmoid2(sum);
            var y2 = -9.8 * d_sigmoid(-4.9 * sum) ;
            console.log(y1, y2);
        </script>
    </body>
</html>